#!/bin/sh
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

# The Step-8 campaign layer is dependency-free.  Run with bytecode disabled so
# this gate does not dirty the source tree on shared validation systems.
PYTHONDONTWRITEBYTECODE=1 python3 "$here/test_campaign.py"
PYTHONDONTWRITEBYTECODE=1 python3 "$here/test_step8_source_contract.py"
