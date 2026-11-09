#!/usr/bin/env python3
"""Check the output of ``amps --validate-brightness`` runs.

Usage (each DIR is the working directory of one fixture run, which contains
``PT/plots``):

  test_i3_brightness_validation.py equivalence DIR_NP1 DIR_NP4 [DIR_NP12 ...]
  test_i3_brightness_validation.py exclude BASE_DIR EXCLUDE_DIR

``equivalence``
  * out=0 (empty sampling window) is identically zero in every column;
  * out=1 values are finite, every impact parameter matches the fixture grid,
    and each column density matches the independent PIC-free chord sum in
    ``brightness_validation_expected_columns.csv``;
  * every numeric column of out=1 agrees between all rank counts.

``exclude``
  BASE_DIR is a run without the out-of-domain particle and EXCLUDE_DIR the
  same decomposition with it under EXCLUDE_AND_REPORT.  Then
  * species without the injected particle are bit-for-bit identical;
  * the H column grows by exactly the injected particle's column (checked via
    the expected-column CSV of each run) and only on rays crossing its shell;
  * excluded column <= added column, excluded counts are positive somewhere,
    and they are zero on rays the particle does not touch;
  * H brightness is unchanged on untouched rays.

The tolerances are relative and only cover floating-point summation order.
"""

import csv
import glob
import math
import os
import re
import sys

REL_TOL = 1.0e-10
TINY = 1.0e-300


def read_tecplot(path):
    names, rows = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith("VARIABLES"):
                names = re.findall(r'"([^"]*)"', line)
            elif line.startswith(("TITLE", "DATASETAUXDATA", "ZONE")):
                continue
            elif line.strip():
                rows.append([float(value) for value in line.split()])
    if names is None or not rows:
        raise AssertionError(f"{path}: missing VARIABLES or data")
    for row in rows:
        if len(row) != len(names):
            raise AssertionError(f"{path}: row has {len(row)} values, expected {len(names)}")
        if not all(math.isfinite(value) for value in row):
            raise AssertionError(f"{path}: non-finite value")
    return names, rows


def brightness_files(run_dir, output_index):
    pattern = os.path.join(run_dir, "PT", "plots", f"pic.*.brightness.s=*.out={output_index}.dat")
    files = sorted(glob.glob(pattern))
    if not files:
        raise AssertionError(f"no brightness tables match {pattern}")
    result = {}
    for path in files:
        match = re.search(r"pic\.(.+)\.brightness\.s=(\d+)\.out=", os.path.basename(path))
        result[(match.group(1), int(match.group(2)))] = path
    return result


def expected_columns(run_dir):
    path = os.path.join(run_dir, "PT", "plots", "brightness_validation_expected_columns.csv")
    table = {}
    with open(path) as handle:
        for row in csv.DictReader(handle):
            table[(int(row["species_index"]), int(row["impact_index"]))] = (
                float(row["impact_parameter_m"]), float(row["expected_column_m2"]))
    return table


def relative_difference(a, b):
    scale = max(abs(a), abs(b))
    return 0.0 if scale < TINY else abs(a - b) / scale


def column_index(names, prefix):
    matches = [i for i, name in enumerate(names) if name.startswith(prefix)]
    if len(matches) != 1:
        raise AssertionError(f"expected exactly one column starting with {prefix!r}: {names}")
    return matches[0]


def check_single_run(run_dir):
    """Empty window, finiteness, grid, and independent column check."""
    worst = 0.0
    for key, path in brightness_files(run_dir, 0).items():
        names, rows = read_tecplot(path)
        for row in rows:
            if any(value != 0.0 for value in row[1:]):
                raise AssertionError(f"{path}: empty window has a nonzero value")
    expected = expected_columns(run_dir)
    for (symbol, species), path in brightness_files(run_dir, 1).items():
        names, rows = read_tecplot(path)
        column = column_index(names, "Column Density")
        nonzero = 0
        for impact, row in enumerate(rows):
            b_expected, n_expected = expected[(species, impact)]
            if relative_difference(row[0], b_expected) > 1.0e-14:
                raise AssertionError(f"{path}: impact parameter {impact} differs from fixture grid")
            diff = relative_difference(row[column], n_expected)
            worst = max(worst, diff)
            if diff > REL_TOL:
                raise AssertionError(
                    f"{path}: ray {impact} column {row[column]:.17e} vs independent {n_expected:.17e}")
            nonzero += row[column] > 0.0
        if nonzero < len(rows) // 2:
            raise AssertionError(f"{path}: only {nonzero}/{len(rows)} rays have a positive column")
        print(f"  {run_dir}: {symbol} rays={len(rows)} positive-columns={nonzero}")
    return worst


def equivalence(run_dirs):
    worst_independent = max(check_single_run(run_dir) for run_dir in run_dirs)
    reference = brightness_files(run_dirs[0], 1)
    worst_rank = 0.0
    for run_dir in run_dirs[1:]:
        candidate = brightness_files(run_dir, 1)
        if set(candidate) != set(reference):
            raise AssertionError(f"{run_dir}: species set differs from {run_dirs[0]}")
        for key, path in reference.items():
            names_a, rows_a = read_tecplot(path)
            names_b, rows_b = read_tecplot(candidate[key])
            if names_a != names_b or len(rows_a) != len(rows_b):
                raise AssertionError(f"{candidate[key]}: layout differs from {path}")
            for row_a, row_b in zip(rows_a, rows_b):
                for a, b in zip(row_a, row_b):
                    diff = relative_difference(a, b)
                    worst_rank = max(worst_rank, diff)
                    if diff > REL_TOL:
                        raise AssertionError(
                            f"{candidate[key]} vs {path}: {a:.17e} != {b:.17e}")
    print(f"PASS brightness decomposition equivalence runs={len(run_dirs)} "
          f"max-rank-rel-diff={worst_rank:.3e} max-independent-column-rel-diff={worst_independent:.3e}")


def exclude(base_dir, exclude_dir):
    base = brightness_files(base_dir, 1)
    excl = brightness_files(exclude_dir, 1)
    expected_base = expected_columns(base_dir)
    expected_excl = expected_columns(exclude_dir)
    if set(base) != set(excl):
        raise AssertionError("species sets differ between base and exclude runs")
    touched_species = 0
    total_excluded = 0
    for key in sorted(base):
        symbol, species = key
        names, rows_base = read_tecplot(base[key])
        names_x, rows_excl = read_tecplot(excl[key])
        col = column_index(names, "Column Density")
        rayleigh = [i for i, n in enumerate(names)
                    if n.startswith("Brightness") and n.endswith("[R]")]
        if len(rayleigh) != 1:
            raise AssertionError(f"{symbol}: expected one Brightness [R] column")
        bright = rayleigh[0]
        frac = column_index(names, "Excluded Column Fraction")
        count = column_index(names, "Excluded g Evaluations")
        added_any = False
        for impact, (row_b, row_x) in enumerate(zip(rows_base, rows_excl)):
            added_expected = (expected_excl[(species, impact)][1] -
                              expected_base[(species, impact)][1])
            added = row_x[col] - row_b[col]
            scale = max(abs(row_x[col]), TINY)
            if abs(added - added_expected) > REL_TOL * scale + 1.0e-300:
                raise AssertionError(f"{symbol} ray {impact}: added column {added:.6e} vs expected {added_expected:.6e}")
            excluded_column = row_x[frac] * row_x[col]
            touched = added_expected > REL_TOL * scale
            if row_b[count] != 0:
                raise AssertionError(f"{symbol} ray {impact}: base run reports excluded evaluations")
            if not touched:
                if row_x[count] != 0:
                    raise AssertionError(f"{symbol} ray {impact}: exclusions on an untouched ray")
                if relative_difference(row_x[bright], row_b[bright]) > REL_TOL:
                    raise AssertionError(f"{symbol} ray {impact}: brightness changed on an untouched ray")
            else:
                added_any = True
                if excluded_column > added * (1.0 + 1.0e-9) + 1.0e-300:
                    raise AssertionError(f"{symbol} ray {impact}: excluded column exceeds added column")
            total_excluded += int(row_x[count])
        if added_any:
            touched_species += 1
        else:
            for row_b, row_x in zip(rows_base, rows_excl):
                if row_b != row_x:
                    raise AssertionError(f"{symbol}: untouched species changed")
    if touched_species != 1:
        raise AssertionError(f"expected exactly one touched species, found {touched_species}")
    if total_excluded <= 0:
        raise AssertionError("EXCLUDE_AND_REPORT recorded no excluded evaluations")
    print(f"PASS brightness EXCLUDE_AND_REPORT accounting excluded-evaluations={total_excluded}")


def main(argv):
    if len(argv) >= 3 and argv[1] == "equivalence":
        equivalence(argv[2:])
    elif len(argv) == 4 and argv[1] == "exclude":
        exclude(argv[2], argv[3])
    else:
        print(__doc__)
        return 2
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv))
    except AssertionError as error:
        print(f"FAIL {error}")
        sys.exit(1)
