"""Deterministic JSON, covariance algebra, and strict source-asset provenance.

Only the standard library is used so preprocessing and its release tests work
on a login node without a scientific Python environment. Covariances are full
matrices, not silently reduced to independent error bars. Publication uses a
sibling temporary directory and one atomic rename, like stage-10 bundles.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import tempfile
from datetime import datetime

class PreprocessingError(ValueError):
    """A typed input/provenance failure; no fallback or partially saved asset."""

def require(condition, message):
    if not condition:
        raise PreprocessingError(message)

def canonical(value):
    """Stable bytes bind identity to content, version and all declared metadata."""
    try:
        return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    except (ValueError, TypeError) as error:
        raise PreprocessingError("noncanonical/nonfinite JSON: " + str(error)) from error

def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()

def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()

def load_json(path):
    # Duplicate keys and JSON's nonstandard NaN/Infinity spellings must fail:
    # accepting them would let two parsers disagree about the frozen authority.
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate JSON key: " + key)
            result[key] = value
        return result
    def constant(value):
        raise PreprocessingError("nonfinite JSON token: " + value)
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"),
                          object_pairs_hook=pairs, parse_constant=constant)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise PreprocessingError("cannot read JSON asset: " + str(error)) from error

def finite(values):
    require(all(isinstance(x, (int, float)) and not isinstance(x, bool)
                and math.isfinite(x) for x in values), "nonfinite/non-numeric values")

def transpose(a):
    require(bool(a) and all(len(row) == len(a[0]) for row in a), "ragged matrix")
    return [list(row) for row in zip(*a)]

def matmul(a, b):
    bt = transpose(b)
    require(all(len(row) == len(b) for row in a), "incompatible matrix dimensions")
    return [[math.fsum(x*y for x, y in zip(row, column)) for column in bt] for row in a]

def matvec(a, x):
    require(all(len(row) == len(x) for row in a), "incompatible vector dimension")
    return [math.fsum(v*w for v, w in zip(row, x)) for row in a]

def covariance(c, size):
    require(len(c) == size and all(len(row) == size for row in c), "covariance dimension mismatch")
    finite(v for row in c for v in row)
    require(all(c[i][i] > 0 for i in range(size)), "covariance has nonpositive variance")
    for i in range(size):
        for j in range(size):
            require(abs(c[i][j]-c[j][i]) <= 1e-12*max(math.sqrt(c[i][i]*c[j][j]), 1e-300),
                    "covariance is not symmetric")
    # Correlation scaling avoids losing an SPD covariance simply because a
    # joint observation mixes meters and dimensionless log-density.
    scale = [math.sqrt(c[i][i]) for i in range(size)]
    correlation = [[c[i][j]/(scale[i]*scale[j]) for j in range(size)] for i in range(size)]
    cholesky(correlation)
    return c

def cholesky(a):
    n = len(a)
    require(n > 0 and all(len(row) == n for row in a), "square nonempty matrix required")
    lower = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            value = a[i][j]-math.fsum(lower[i][k]*lower[j][k] for k in range(j))
            if i == j:
                require(value > 1e-14, "covariance/design matrix is not positive definite")
                lower[i][j] = math.sqrt(value)
            else:
                lower[i][j] = value/lower[j][j]
    return lower

def solve_spd(a, b):
    require(len(a) == len(b) and bool(a) and all(len(row) == len(a) for row in a), "linear system dimension mismatch")
    require(all(a[i][i] > 0 and math.isfinite(a[i][i]) for i in range(len(a))), "singular/nonpositive design diagonal")
    scale = [math.sqrt(a[i][i]) for i in range(len(a))]
    lower = cholesky([[a[i][j]/(scale[i]*scale[j]) for j in range(len(a))] for i in range(len(a))])
    y = []
    for i in range(len(a)):
        y.append((b[i]/scale[i]-math.fsum(lower[i][j]*y[j] for j in range(i)))/lower[i][i])
    x = [0.0]*len(a)
    for i in range(len(a)-1, -1, -1):
        x[i] = (y[i]-math.fsum(lower[j][i]*x[j] for j in range(i+1, len(a))))/lower[i][i]
    return [x[i]/scale[i] for i in range(len(a))]

def inverse_spd(a):
    n = len(a)
    result = transpose([solve_spd(a, [float(i == j) for i in range(n)]) for j in range(n)])
    return [[0.5*(result[i][j]+result[j][i]) for j in range(n)] for i in range(n)]

def gaussian(residual, c):
    """Full normalized Gaussian likelihood, including covariance determinant."""
    finite(residual)
    covariance(c, len(residual))
    scale = [math.sqrt(c[i][i]) for i in range(len(c))]
    lower = cholesky([[c[i][j]/(scale[i]*scale[j]) for j in range(len(c))] for i in range(len(c))])
    logdet = 2*math.fsum(math.log(x) for x in scale)+2*math.fsum(math.log(lower[i][i]) for i in range(len(c)))
    chi2 = math.fsum(x*y for x, y in zip(residual, solve_spd(c, residual)))
    return {"chi2": chi2, "log_likelihood": -0.5*(chi2+logdet+len(c)*math.log(2*math.pi))}

def linear_fit(design, values, c):
    """Generalized least squares; parameter covariance is (A^T C^-1 A)^-1."""
    finite(values)
    covariance(c, len(values))
    require(len(design) == len(values), "fit sample dimension mismatch")
    finite(v for row in design for v in row)
    at = transpose(design)
    weighted_columns = transpose([solve_spd(c, column) for column in at])
    normal = matmul(at, weighted_columns)
    parameter_covariance = inverse_spd(normal)
    parameters = solve_spd(normal, matvec(at, solve_spd(c, values)))
    predicted = matvec(design, parameters)
    residual = [y-p for y, p in zip(values, predicted)]
    return {"parameters": parameters, "covariance": parameter_covariance,
            "chi2": gaussian(residual, c)["chi2"], "residual": residual}

def propagate(jacobian, c):
    result = matmul(matmul(jacobian, c), transpose(jacobian))
    # Mathematical covariance is symmetric; compute/retain that symmetry
    # explicitly rather than exporting roundoff from two differently ordered
    # matrix dot products, especially for ill-scaled geometric parameters.
    return [[0.5*(result[i][j]+result[j][i]) for j in range(len(result))] for i in range(len(result))]

ROLES = {"construction", "qualification", "withheld-validation"}
METADATA = {"asset_id", "kind", "units", "frame", "epoch_utc", "processing_version",
            "data_use_role", "source_uri", "source_sha256", "coordinate_definition",
            "support", "independence_id"}

def validate_metadata(m):
    require(METADATA <= set(m), "missing asset metadata: " + ", ".join(sorted(METADATA-set(m))))
    for key in METADATA-{"support"}:
        require(isinstance(m[key], str) and bool(m[key]), "empty/non-string metadata: " + key)
    require(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", m["asset_id"]) is not None, "unsafe asset ID")
    require(re.fullmatch(r"[a-f0-9]{64}", m["source_sha256"]) is not None, "invalid source checksum")
    require(m["data_use_role"] in ROLES, "unsupported data-use role")
    require(re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", m["epoch_utc"]) is not None, "epoch must be explicit UTC ISO8601")
    try:
        datetime.strptime(m["epoch_utc"], "%Y-%m-%dT%H:%M:%SZ")
    except ValueError as error:
        raise PreprocessingError("invalid UTC calendar epoch") from error
    require(m["units"] != "unknown" and m["frame"] != "unknown", "ambiguous units/frame")
    require(isinstance(m["support"], dict) and bool(m["support"]), "missing declared sample support")
    canonical(m)

def validate_asset_roles(assets):
    require(bool(assets), "empty observation asset bundle")
    ids = set()
    authorities = {}
    epochs = set()
    versions = set()
    for asset in assets:
        m = asset["metadata"]
        validate_metadata(m)
        require(m["asset_id"] not in ids, "duplicate asset ID")
        ids.add(m["asset_id"]); epochs.add(m["epoch_utc"]); versions.add(m["processing_version"])
        # Renaming an observation, response, or source path cannot manufacture
        # independent evidence. Reuse across roles is detected by BYTES and by
        # the declared independence/response groups, not filename alone.
        for category, value in [("bytes", m["source_sha256"]), ("group", m["independence_id"]),
                                ("response", m.get("response_id"))]:
            if value:
                key = (category, value)
                role = m["data_use_role"]
                require(key not in authorities or authorities[key] == role,
                        "construction/qualification/withheld response or observation reuse")
                authorities[key] = role
    require(len(epochs) == 1, "mixed reference epochs in asset bundle")
    require(len(versions) == 1, "mixed processing versions in asset bundle")

def freeze_record(record):
    require("identity" not in record, "record is already frozen")
    result = dict(record)
    result["identity"] = digest(record)
    return result

def verify_frozen(record):
    require(isinstance(record, dict) and "identity" in record, "missing immutable identity")
    plain = {k: v for k, v in record.items() if k != "identity"}
    require(record["identity"] == digest(plain), "post hoc change to frozen record")
    return record

def publish_bundle(destination, assets, sources, processing_version):
    validate_asset_roles(assets)
    require(all(a["metadata"]["source_sha256"] in sources for a in assets), "missing source bytes at publication")
    require(all(a["metadata"]["processing_version"] == processing_version for a in assets), "processing version differs from bundle")
    target = Path(destination)
    require(not target.exists(), "immutable output already exists")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix="."+target.name+"-", dir=str(target.parent)))
    try:
        members = []
        (temporary/"assets").mkdir(); (temporary/"sources").mkdir()
        for asset in sorted(assets, key=lambda a: a["metadata"]["asset_id"]):
            data = canonical(asset)
            name = "assets/"+asset["metadata"]["asset_id"]+".json"
            (temporary/name).write_bytes(data)
            members.append({"path": name, "sha256": sha256_bytes(data)})
        for checksum, source in sorted(sources.items()):
            require(sha256_bytes(source) == checksum, "source bytes differ from metadata")
            name = "sources/"+checksum+".json"
            (temporary/name).write_bytes(source)
            members.append({"path": name, "sha256": checksum})
        manifest = freeze_record({"schema": "sep-observation-assets-v1",
            "processing_version": processing_version,
            "epoch_utc": assets[0]["metadata"]["epoch_utc"], "members": members})
        (temporary/"manifest.json").write_bytes(canonical(manifest))
        temporary.rename(target)
        return manifest
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)

def read_bundle(directory):
    base = Path(directory)
    manifest = verify_frozen(load_json(base/"manifest.json"))
    require(manifest["schema"] == "sep-observation-assets-v1", "unsupported observation bundle schema")
    assets = []
    names = set()
    require(manifest["members"], "empty observation manifest")
    for member in manifest["members"]:
        name = member["path"]
        require(isinstance(name, str) and not Path(name).is_absolute() and
                ".." not in Path(name).parts and name not in names, "unsafe/duplicate manifest member")
        names.add(name)
        path = base/name
        require(path.is_file() and not path.is_symlink() and base.resolve() in path.resolve().parents, "missing/symlink/outside asset member")
        require(sha256_bytes(path.read_bytes()) == member["sha256"], "asset/source member checksum changed")
        if name.startswith("assets/"):
            assets.append(load_json(path))
    validate_asset_roles(assets)
    require(all(a["metadata"]["epoch_utc"] == manifest["epoch_utc"] and
                a["metadata"]["processing_version"] == manifest["processing_version"] for a in assets),
            "bundle/asset epoch or processing metadata differ")
    for asset in assets:
        require("sources/"+asset["metadata"]["source_sha256"]+".json" in names,
                "source bytes absent from immutable bundle")
    return manifest, assets
