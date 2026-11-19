"""Offline query CLI; every command uses an explicitly pinned snapshot.

Human-readable help documents that source integrity and physical validation
are different operations. Output is JSON unless --csv requests exact table
cells. No CLI operation rewrites source data or chooses a physical model.
"""

import argparse
import csv
import json
from pathlib import Path
import sys

from .api import Database, DatabaseError, plain


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, required=True,
                        help="database-lock.json for an exact installed release/profile/export")
    commands = parser.add_subparsers(dest="command", required=True)
    listing = commands.add_parser("list", help="list source datasets, tables or source objects")
    listing.add_argument("--kind", choices=("datasets", "tables", "sources", "species", "reactions", "models"), default="datasets")
    query = commands.add_parser("query", help="search dataset metadata without source ranking")
    query.add_argument("text")
    describe = commands.add_parser("describe", help="show dataset/table/source provenance and limitations")
    describe.add_argument("identifier")
    table = commands.add_parser("table", help="read a permission-cleared CSV; preserve native cells/order")
    table.add_argument("identifier")
    table.add_argument("--csv", action="store_true", help="emit native table cells instead of JSON")
    source = commands.add_parser("source-path", help="resolve a cleared immutable source object")
    source.add_argument("identifier")
    verify = commands.add_parser("verify", help="verify archive bytes/schema; NOT physical validation")
    verify.add_argument("--shallow", action="store_true")
    commands.add_parser("identity", help="print release/profile/runtime-export pins")
    args = parser.parse_args(argv)
    try:
        database = Database.from_lock(args.lock)
        if args.command == "list":
            result = getattr(database, args.kind)()
        elif args.command == "query":
            result = database.query(args.text)
        elif args.command == "describe":
            result = database.describe(args.identifier)
        elif args.command == "source-path":
            result = {"source_id": args.identifier, "path": str(database.source_path(args.identifier))}
        elif args.command == "verify":
            result = database.verify(deep=not args.shallow)
        elif args.command == "identity":
            result = database.identity
        else:
            value = database.table(args.identifier)
            if args.csv:
                writer = csv.writer(sys.stdout, lineterminator="\n")
                writer.writerow(value.columns)
                writer.writerows(value.rows)
                return 0
            result = {"columns": value.columns, "rows": value.rows,
                      "metadata": value.metadata, "row_count": value.row_count}
        print(json.dumps(plain(result), indent=2, ensure_ascii=False, allow_nan=False))
        return 0
    except DatabaseError as error:
        print(json.dumps({"status": "SKIPPED" if error.code in ("Unavailable", "PermissionRequired") else "ERROR",
                          "code": error.code, "message": str(error)}, indent=2), file=sys.stderr)
        return 77 if error.code in ("Unavailable", "PermissionRequired") else 2
    except (OSError, ValueError, KeyError) as error:
        print(json.dumps({"status": "ERROR", "code": "InfrastructureError", "message": str(error)}, indent=2), file=sys.stderr)
        return 2
