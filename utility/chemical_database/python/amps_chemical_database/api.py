"""Strict, source-preserving Python API for an immutable chemical release.

Units and frames are metadata attached to particular source assertions, never
guessed from a species name or a filename. CSV cells remain original strings.
The reader checks the same runtime source/table indexes as the C++ Catalogue;
no server, network connection or implicit latest-release selection is used.
"""

import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json
from pathlib import Path
import re
import sys
from types import MappingProxyType


SOURCE_COLUMNS = (
    "source_id", "path", "sha256", "media_type", "availability",
    "storage_permission", "evaluation_permission", "capabilities",
    "source_uri", "notes",
)
TABLE_COLUMNS = (
    "table_id", "path", "sha256", "format", "schema_id", "source_ids",
    "availability", "evaluation_permission", "capabilities", "columns",
    "row_count", "notes",
)
EXPORT_KEYS = frozenset((
    "schema", "codec", "release_id", "release_manifest_path",
    "release_manifest_sha256", "profile_id", "profile_path",
    "profile_sha256", "source_index_path", "source_index_sha256",
    "table_index_path", "table_index_sha256", "profile_capabilities",
    "topology_fingerprint",
))
PERMISSIONS = frozenset(("allowed", "pending", "prohibited", "unknown", "not-applicable"))
AVAILABILITIES = frozenset(("available", "unavailable", "metadata-only", "quarantined"))
_DIGEST = re.compile(r"[0-9a-f]{64}\Z")
_IDENTIFIER = re.compile(r"[A-Za-z0-9_.:-]+\Z")


class DatabaseError(RuntimeError):
    """Typed failure; missing data never become a numerical zero.

    ``code`` matches the shared C++ outcome vocabulary. Test runners should
    distinguish unavailable prerequisites from corrupt input/infrastructure.
    Production callers decide whether a required missing capability is fatal.
    """

    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


def _fail(code, message):
    raise DatabaseError(code, message)


def sha256(path):
    """Hash a file in bounded I/O blocks; do not load large native objects."""
    digest = hashlib.sha256()
    try:
        with Path(path).open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        _fail("Unavailable", "cannot read source object {}: {}".format(path, error))
    return digest.hexdigest()


def _digest(value, label):
    if not isinstance(value, str) or not _DIGEST.fullmatch(value):
        _fail("InvalidInput", "{} must be a lowercase SHA-256 digest".format(label))
    return value


def member(root, relative):
    """Resolve only root-contained paths, including symlink targets.

    Cross-machine exports use POSIX relative paths. Absolute paths, parent
    components and backslashes are refused before accessing the filesystem.
    Historical external original paths may exist in documentary JSON, but
    the runtime reader never follows them to find a missing imported object.
    """
    if not isinstance(relative, str) or not relative or "\\" in relative or "\x00" in relative:
        _fail("CorruptSource", "empty or nonportable catalogue member path")
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts:
        _fail("CorruptSource", "member path escapes database root: " + relative)
    resolved = (Path(root) / path).resolve()
    try:
        resolved.relative_to(Path(root).resolve())
    except ValueError:
        _fail("CorruptSource", "member symlink escapes database root: " + relative)
    return resolved


def strict_json(path):
    """Reject duplicate fields and nonstandard nonfinite JSON constants."""
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                _fail("CorruptSource", "duplicate JSON key {} in {}".format(key, path))
            result[key] = value
        return result

    def constant(value):
        _fail("CorruptSource", "nonfinite JSON token {} in {}".format(value, path))

    try:
        return json.loads(Path(path).read_text(encoding="utf-8"),
                          object_pairs_hook=pairs, parse_constant=constant)
    except (ValueError, UnicodeError, OSError) as error:
        _fail("CorruptSource", "invalid JSON {}: {}".format(path, error))


def _mapping(value, label):
    if not isinstance(value, dict):
        _fail("CorruptSource", label + " must be a JSON object")
    return value


def _freeze(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def plain(value):
    """Convert defensive API views into JSON-serializable plain objects."""
    if isinstance(value, (dict, MappingProxyType)):
        return {key: plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(item) for item in value]
    return value


def _split(value):
    if value in ("", "none"):
        return ()
    parts = tuple(value.split(";"))
    if any(not _IDENTIFIER.fullmatch(item) for item in parts) or len(set(parts)) != len(parts):
        _fail("CorruptSource", "invalid or duplicate semicolon-list entry: " + value)
    return parts


def _columns(value):
    """Decode exact native header labels, independently of identifier lists.

    Native headers may contain significant padding, spaces and punctuation.
    Capability IDs forbid padding; applying that policy to column labels would
    silently change a scientific source schema or reject a valid source table.
    """
    parts = tuple(value.split(";")) if value else ()
    if not parts or any(not item or "\r" in item or "\n" in item for item in parts) or len(set(parts)) != len(parts):
        _fail("CorruptSource", "invalid exact CSV column list")
    return parts


def _verified(root, relative, expected):
    _digest(expected, relative)
    path = member(root, relative)
    if not path.is_file():
        _fail("Unavailable", "required database object is absent: " + relative)
    actual = sha256(path)
    if actual != expected:
        _fail("CorruptSource", "SHA-256 mismatch for {}: expected {}, actual {}".format(
            relative, expected, actual))
    return path


def _read_csv(path, expected_columns=None):
    try:
        with Path(path).open("r", encoding="utf-8", newline="") as stream:
            text = stream.read()
        if "\x00" in text:
            _fail("CorruptSource", "CSV contains a NUL byte: " + str(path))
        # Python's csv dialect accepts a quote embedded in an unquoted field.
        # The scientific archive does not: check the shared strict grammar
        # before parsing so Python and native consumers reject the same bytes.
        quoted = closed = started = False
        i = 0
        while i < len(text):
            char = text[i]
            if quoted:
                if char == '"':
                    if i + 1 < len(text) and text[i + 1] == '"':
                        i += 1
                    else:
                        quoted, closed = False, True
            elif char in (",", "\n", "\r"):
                if char == "\r":
                    if i + 1 == len(text) or text[i + 1] != "\n":
                        _fail("CorruptSource", "bare CR outside a quoted CSV field")
                    i += 1
                closed = started = False
            elif char == '"':
                if started or closed:
                    _fail("CorruptSource", "quote inside an unquoted CSV field")
                quoted, started = True, True
            else:
                if closed:
                    _fail("CorruptSource", "characters after a closing CSV quote")
                started = True
            i += 1
        if quoted:
            _fail("CorruptSource", "unterminated CSV quote")
        # The C++ reader implements the same quote/field/row rules. The strict
        # dialect matters: malformed quoting must not alter source references.
        rows = list(csv.reader(io.StringIO(text, newline=""), strict=True))
    except (csv.Error, UnicodeError, OSError) as error:
        _fail("CorruptSource", "cannot parse CSV {}: {}".format(path, error))
    if not rows or not rows[0] or any(not field for field in rows[0]):
        _fail("CorruptSource", "CSV has no complete header: " + str(path))
    columns = tuple(rows[0])
    if len(set(columns)) != len(columns):
        _fail("CorruptSource", "duplicate CSV header: " + str(path))
    if expected_columns is not None and columns != tuple(expected_columns):
        _fail("CorruptSource", "CSV header/schema mismatch: " + str(path))
    if any(len(row) != len(columns) for row in rows[1:]):
        _fail("CorruptSource", "CSV row has an inconsistent field count: " + str(path))
    return columns, tuple(tuple(row) for row in rows[1:])


def _export(path):
    """Parse the fixed, flat projection of the existing AMPS block grammar.

    Runtime exports are produced by the shared publisher, not general user
    input. They intentionally forbid includes, continuations and child blocks.
    The general C++ grammar reader applies the matching owner restrictions.
    """
    try:
        lines = Path(path).read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        _fail("CorruptSource", "cannot read runtime export: " + str(error))
    result = {}
    opened = closed = False
    for number, physical in enumerate(lines, 1):
        line = physical.split("!", 1)[0].strip()
        if not line:
            continue
        if line == "#begin chemical_database: runtime" and not opened and not closed:
            opened = True
            continue
        if line == "#end chemical_database" and opened and not closed:
            closed = True
            continue
        if not opened or closed or line.startswith("#") or "=" not in line or line.endswith("\\"):
            _fail("CorruptSource", "invalid export structure at line {}".format(number))
        key, value = (field.strip() for field in line.split("=", 1))
        if key not in EXPORT_KEYS or key in result or not value:
            _fail("CorruptSource", "unknown/duplicate/empty export field at line {}".format(number))
        result[key] = value
    if not closed or set(result) != EXPORT_KEYS:
        _fail("CorruptSource", "incomplete runtime export")
    if result["schema"] != "amps-chemical-runtime-v1" or result["codec"] != "csv-block-v1":
        _fail("UnsupportedSchema", "unsupported runtime export schema/codec")
    if result["topology_fingerprint"] != "none":
        _digest(result["topology_fingerprint"], "topology_fingerprint")
    _split(result["profile_capabilities"])
    return result


@dataclass(frozen=True)
class Table:
    """Owned lexical CSV table; row order and missing cells are preserved."""
    columns: tuple
    rows: tuple
    metadata: object

    @property
    def row_count(self):
        return len(self.rows)

    def dictionaries(self):
        for row in self.rows:
            yield MappingProxyType(dict(zip(self.columns, row)))

    def value(self, row, column):
        try:
            index = self.columns.index(column)
        except ValueError:
            _fail("InvalidInput", "unknown table column: " + column)
        return self.rows[row][index]

    def number(self, row, column, expected_unit):
        """Return an exact Decimal only after an explicit unit agreement.

        There is no implicit unit conversion or parsing of native qualifiers,
        uncertainty parentheses or fill values as physical numerical values.
        """
        units = self.metadata.get("units", {})
        if not expected_unit or units.get(column) != expected_unit:
            _fail("Unavailable", "column unit is absent or differs from requested unit")
        text = self.value(row, column)
        if not text:
            _fail("Unavailable", "source numerical cell is missing")
        try:
            value = Decimal(text)
        except InvalidOperation:
            _fail("Unavailable", "cell contains a source qualifier or unsupported numerical form")
        if not value.is_finite():
            _fail("CorruptSource", "nonfinite numerical source cell")
        return value


class Database:
    """A hash-pinned immutable catalogue and its runtime selection.

    ``from_lock`` is convenient for the explicitly installed snapshot. A lock
    identifies an exact release/profile/export; it never means 'latest'.
    Scientific use still requires the relevant model/application qualification.
    """

    def __init__(self):
        _fail("InvalidInput", "use Database.open or Database.from_lock")

    @classmethod
    def from_lock(cls, lock_path):
        lock_path = Path(lock_path).resolve()
        lock = _mapping(strict_json(lock_path), "installation lock")
        if lock.get("schema") != "amps-chemical-database-lock-v1":
            _fail("UnsupportedSchema", "unknown installation lock schema")
        required = ("release_id", "release_manifest_sha256", "profile_id",
                    "profile_sha256", "export_path", "export_manifest_sha256")
        if any(not isinstance(lock.get(key), str) or not lock[key] for key in required):
            _fail("CorruptSource", "installation lock is missing a required identity field")
        # Installation locks are root-local and portable. They may not redirect
        # the reader to another private/home source directory.
        if lock.get("database_root", ".") != ".":
            _fail("InvalidInput", "installation lock database_root must be '.'")
        return cls.open(lock_path.parent, lock["release_id"], lock["release_manifest_sha256"],
                        lock["profile_id"], lock["profile_sha256"],
                        lock["export_path"], lock["export_manifest_sha256"])

    @classmethod
    def open(cls, root, release_id, expected_manifest_sha256, profile_id,
             expected_profile_sha256, export_path, expected_export_sha256):
        root = Path(root).resolve()
        export = _export(_verified(root, export_path, expected_export_sha256))
        for key, expected in (("release_id", release_id),
                              ("release_manifest_sha256", expected_manifest_sha256),
                              ("profile_id", profile_id),
                              ("profile_sha256", expected_profile_sha256)):
            if export[key] != expected:
                _fail("CorruptSource", "runtime export disagrees with pinned " + key)
        release_path = _verified(root, export["release_manifest_path"], expected_manifest_sha256)
        profile_path = _verified(root, export["profile_path"], expected_profile_sha256)
        release = _mapping(strict_json(release_path), "release manifest")
        profile = _mapping(strict_json(profile_path), "profile")
        if release.get("schema") != "amps-chemical-release-v1" or release.get("release_id") != release_id:
            _fail("UnsupportedSchema", "canonical release schema/identity mismatch")
        if profile.get("schema") != "amps-chemical-profile-v1" or profile.get("id") != profile_id:
            _fail("UnsupportedSchema", "canonical profile schema/identity mismatch")
        # Profiles are members of the pinned canonical release. An attacker
        # cannot substitute another JSON file simply by updating the export.
        members = _mapping(release.get("members"), "release members")
        if members.get(export["profile_path"]) != expected_profile_sha256:
            _fail("CorruptSource", "profile is not a member of the canonical release")
        sources_path = _verified(root, export["source_index_path"], export["source_index_sha256"])
        tables_path = _verified(root, export["table_index_path"], export["table_index_sha256"])
        _, source_rows = _read_csv(sources_path, SOURCE_COLUMNS)
        _, table_rows = _read_csv(tables_path, TABLE_COLUMNS)
        sources, tables = {}, {}
        for row in source_rows:
            source = dict(zip(SOURCE_COLUMNS, row))
            key = source["source_id"]
            if not _IDENTIFIER.fullmatch(key) or key in sources:
                _fail("CorruptSource", "duplicate or empty source identifier")
            cls._validate_source(source, root)
            sources[key] = source
        for row in table_rows:
            table = dict(zip(TABLE_COLUMNS, row))
            key = table["table_id"]
            if not _IDENTIFIER.fullmatch(key) or key in tables:
                _fail("CorruptSource", "duplicate or empty table identifier")
            cls._validate_table(table, sources, root)
            tables[key] = table
        catalogue_ref = _mapping(release.get("catalogue"), "release catalogue reference")
        catalogue_path = _verified(root, catalogue_ref.get("path", ""), catalogue_ref.get("sha256", ""))
        if members.get(catalogue_ref["path"]) != catalogue_ref["sha256"]:
            _fail("CorruptSource", "catalogue is not a release member")
        catalogue = _mapping(strict_json(catalogue_path), "catalogue")
        if catalogue.get("schema") != "amps-chemical-catalogue-v1":
            _fail("UnsupportedSchema", "unknown canonical catalogue schema")
        datasets = catalogue.get("datasets", [])
        if not isinstance(datasets, list) or any(not isinstance(item, dict) or
                not isinstance(item.get("id"), str) or not item["id"] for item in datasets):
            _fail("CorruptSource", "catalogue datasets must be named JSON objects")
        _mapping(catalogue.get("table_metadata", {}), "table metadata")
        if len({item["id"] for item in datasets}) != len(datasets):
            _fail("CorruptSource", "duplicate canonical dataset identifier")
        # Construct only after every check; a partially loaded object is never
        # published. Stored metadata is frozen recursively and queries return
        # defensive immutable views.
        self = object.__new__(cls)
        self.root = root
        self.identity = _freeze(dict(export, export_path=export_path,
                                     export_sha256=expected_export_sha256))
        self._release = _freeze(release)
        self._profile = _freeze(profile)
        self._catalogue = _freeze(catalogue)
        self._sources = _freeze(sources)
        self._tables = _freeze(tables)
        self._datasets = _freeze({item["id"]: item for item in datasets})
        return self

    @staticmethod
    def _validate_source(source, root):
        if source["availability"] not in AVAILABILITIES:
            _fail("CorruptSource", "unknown source availability")
        for field in ("storage_permission", "evaluation_permission"):
            if source[field] not in PERMISSIONS:
                _fail("CorruptSource", "unknown source permission")
        _split(source["capabilities"])
        present = bool(source["path"])
        if present != bool(source["sha256"]):
            _fail("CorruptSource", "source payload path/hash missing its pair")
        if present:
            member(root, source["path"])
            _digest(source["sha256"], source["source_id"])
        if source["availability"] == "available" and not present:
            _fail("CorruptSource", "available source has no payload")
        if present and source["storage_permission"] != "allowed":
            _fail("PermissionRequired", "source payload lacks managed-storage permission")

    @staticmethod
    def _validate_table(table, sources, root):
        if table["availability"] not in AVAILABILITIES or table["evaluation_permission"] not in PERMISSIONS:
            _fail("CorruptSource", "unknown table availability/permission")
        if table["format"] not in ("csv", "native-ascii", "opaque"):
            _fail("UnsupportedSchema", "unknown runtime table format")
        _split(table["capabilities"])
        dependencies = _split(table["source_ids"])
        if not dependencies or any(key not in sources for key in dependencies):
            _fail("CorruptSource", "table refers to absent source identifiers")
        if bool(table["path"]) != bool(table["sha256"]):
            _fail("CorruptSource", "table payload path/hash missing its pair")
        if table["path"]:
            member(root, table["path"])
            _digest(table["sha256"], table["table_id"])
        if table["availability"] == "available" and not table["path"]:
            _fail("CorruptSource", "available table has no payload")
        if table["format"] == "csv":
            _columns(table["columns"])
            if not re.fullmatch(r"[0-9]+", table["row_count"]) or int(table["row_count"]) > sys.maxsize:
                _fail("CorruptSource", "CSV row count must be an exact nonnegative integer")
        elif table["columns"] or table["row_count"] != "unknown":
            _fail("CorruptSource", "opaque/native table has invented CSV shape")

    def datasets(self):
        return tuple(self._datasets.values())

    def species(self):
        return tuple(self._catalogue.get("species", ()))

    def reactions(self):
        return tuple(self._catalogue.get("reactions", ()))

    def sources(self):
        return tuple(self._sources.values())

    def tables(self):
        return tuple(self._table_description(key) for key in self._tables)

    def _table_description(self, key):
        details = plain(self._catalogue.get("table_metadata", {}).get(key, {}))
        if not isinstance(details, dict):
            _fail("CorruptSource", "canonical table details must be an object")
        # Documentary details may add units/name/context but never replace the
        # authenticated runtime path, permission, schema, shape or identifier.
        if any(field in self._tables[key] and details[field] != self._tables[key][field]
               for field in details):
            _fail("CorruptSource", "canonical table details conflict with runtime index")
        return _freeze(dict(plain(self._tables[key]), **details))

    def _resolve_table_id(self, identifier):
        if identifier in self._tables:
            return identifier
        matches = [key for key in self._tables
                   if self._catalogue.get("table_metadata", {}).get(key, {}).get("name") == identifier]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            _fail("InvalidInput", "table name is ambiguous; use its complete stable ID")
        _fail("Unavailable", "table is not in this runtime selection: " + identifier)

    def describe(self, identifier):
        if identifier in self._tables:
            return self._table_description(identifier)
        for mapping in (self._datasets, self._tables, self._sources):
            if identifier in mapping:
                return mapping[identifier]
        return self._table_description(self._resolve_table_id(identifier))

    def query(self, text):
        """Case-insensitive metadata search; never selects a preferred model."""
        needle = text.casefold()
        return tuple(item for item in self.datasets()
                     if needle in json.dumps(plain(item), sort_keys=True).casefold())

    def require(self, capability):
        if not isinstance(capability, str) or not _IDENTIFIER.fullmatch(capability):
            _fail("InvalidInput", "capability identifier must be explicit")
        if capability not in _split(self.identity["profile_capabilities"]):
            _fail("Unavailable", "profile does not declare capability: " + capability)

    def source_path(self, identifier, capability="source-access"):
        try:
            source = self._sources[identifier]
        except KeyError:
            _fail("Unavailable", "source is not in this runtime selection: " + identifier)
        self._admit(source, capability)
        return _verified(self.root, source["path"], source["sha256"])

    def _admit(self, record, capability):
        self.require(capability)
        if capability not in _split(record["capabilities"]):
            _fail("Unavailable", "source/table lacks requested capability: " + capability)
        if "storage_permission" in record and record["storage_permission"] != "allowed":
            _fail("PermissionRequired", "source managed-storage operation is not cleared")
        if record["evaluation_permission"] != "allowed":
            _fail("PermissionRequired", "source/table operation is not cleared: " +
                  record.get("source_id", record.get("table_id", "")))
        if record["availability"] != "available":
            _fail("Unavailable", "source/table is " + record["availability"])

    def table(self, identifier, capability="table-access"):
        identifier = self._resolve_table_id(identifier)
        try:
            record = self._tables[identifier]
        except KeyError:
            _fail("Unavailable", "table is not in this runtime selection: " + identifier)
        self._admit(record, capability)
        for source_id in _split(record["source_ids"]):
            source = self._sources[source_id]
            if source["evaluation_permission"] != "allowed" or source["storage_permission"] != "allowed":
                _fail("PermissionRequired", "table dependency is not cleared: " + source_id)
            if source["availability"] != "available":
                _fail("Unavailable", "table dependency is unavailable: " + source_id)
            _verified(self.root, source["path"], source["sha256"])
        if record["format"] != "csv":
            _fail("UnsupportedSchema", "table requires its native source reader, not the CSV API")
        path = _verified(self.root, record["path"], record["sha256"])
        columns, rows = _read_csv(path, _columns(record["columns"]))
        if len(rows) != int(record["row_count"]):
            _fail("CorruptSource", "runtime table row count changed")
        return Table(columns, rows, self._table_description(identifier))

    def verify(self, deep=True):
        """Verify stored bytes and CSV shape, not physical validation/readiness."""
        checked = 0
        for relative, digest in self._release["members"].items():
            _verified(self.root, relative, digest)
            checked += 1
        if deep:
            for source in self._sources.values():
                if source["path"]:
                    _verified(self.root, source["path"], source["sha256"])
                    checked += 1
            for identifier, table in self._tables.items():
                if table["path"]:
                    path = _verified(self.root, table["path"], table["sha256"])
                    checked += 1
                    if table["format"] == "csv":
                        _, rows = _read_csv(path, _columns(table["columns"]))
                        if len(rows) != int(table["row_count"]):
                            _fail("CorruptSource", "table shape mismatch: " + identifier)
        return {"status": "PASS", "scope": "stored byte/schema integrity only",
                "objects_checked": checked, "datasets": len(self._datasets),
                "sources": len(self._sources), "tables": len(self._tables),
                "available_tables": sum(row["availability"] == "available" for row in self._tables.values()),
                "scientific_validation": "NOT EVALUATED", "physical_activation": False}
