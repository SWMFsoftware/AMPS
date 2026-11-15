"""Test-family registry (utility/nightly/families.json)."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

DEFAULT_PATH = Path(__file__).resolve().parent.parent / "families.json"
UNASSIGNED = "unassigned"
OTHER_GROUP = "other"
_ID = re.compile(r"^[a-z0-9][a-z0-9-]{0,62}$")


class RegistryError(ValueError):
    """families.json is malformed."""


class Registry:
    """Ordered groups and families; the report shows them in this order."""

    def __init__(self, data: Dict[str, Any]):
        if not isinstance(data, dict) or data.get("schema_version") != 1:
            raise RegistryError("families.json must be an object with schema_version 1")
        self.groups: List[Dict[str, str]] = []
        self.families: List[Dict[str, Any]] = []
        group_ids = set()
        for group in data.get("groups", []):
            gid = group.get("id", "")
            if not _ID.match(gid) or gid in group_ids:
                raise RegistryError("bad or duplicate group id: {!r}".format(gid))
            group_ids.add(gid)
            self.groups.append({"id": gid, "title": str(group.get("title", gid)),
                                "source": str(group.get("source", ""))})
        if OTHER_GROUP not in group_ids:
            self.groups.append({"id": OTHER_GROUP, "title": "Other results", "source": ""})
            group_ids.add(OTHER_GROUP)
        seen = set()
        for family in data.get("families", []):
            fid = family.get("id", "")
            if not _ID.match(fid) or fid in seen:
                raise RegistryError("bad or duplicate family id: {!r}".format(fid))
            if family.get("group") not in group_ids:
                raise RegistryError("family {} names unknown group {!r}".format(fid, family.get("group")))
            prefixes = family.get("suite_prefixes", [])
            if not isinstance(prefixes, list) or not all(isinstance(p, str) and p for p in prefixes):
                raise RegistryError("family {}: suite_prefixes must be a list of strings".format(fid))
            seen.add(fid)
            self.families.append({"id": fid, "group": family["group"],
                                  "title": str(family.get("title", fid)),
                                  "runner": str(family.get("runner", "")),
                                  "suite_prefixes": list(prefixes)})
        self._by_id = {f["id"]: f for f in self.families}
        self._groups_by_id = {g["id"]: g for g in self.groups}

    @classmethod
    def load(cls, path: Optional[Path] = None) -> "Registry":
        path = Path(path) if path else DEFAULT_PATH
        try:
            return cls(json.loads(path.read_text()))
        except (OSError, ValueError) as error:
            raise RegistryError("cannot read {}: {}".format(path, error))

    def family(self, fid: str) -> Optional[Dict[str, Any]]:
        return self._by_id.get(fid)

    def group(self, gid: str) -> Dict[str, str]:
        return self._groups_by_id.get(gid) or self._groups_by_id[OTHER_GROUP]

    def family_for_suite(self, suite_id: str) -> str:
        """Family of a src/test suite or application campaign ID; first matching prefix wins."""
        for family in self.families:
            if any(suite_id.startswith(p) for p in family["suite_prefixes"]):
                return family["id"]
        return UNASSIGNED

    def describe(self, fid: str) -> Dict[str, Any]:
        """Family record for display; unregistered IDs are shown in the 'other' group."""
        found = self._by_id.get(fid)
        if found:
            return found
        return {"id": fid, "group": OTHER_GROUP, "title": fid, "runner": "", "suite_prefixes": []}
