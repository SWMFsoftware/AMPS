#!/usr/bin/env python3
"""Build a time-ordered COPS housekeeping table from the ROSINA COPS L2 'bg' files.

DATA SOURCE
  Archive    ROSINA team (PI K. Altwegg, Univ. Bern), PDS3 data set RO-C-ROSINA-2-EXT3-V2.0 ("ROSETTA-ORBITER 67P
             ROSINA 2 EXT3 V2.0"), CODMAC level 2 (raw data as received from the instrument, pressures calibrated by
             the on-board DPU for N2 at 20 C), product creation 2019-02.
             ESA PSA: https://archives.esac.esa.int/psa/ftp/INTERNATIONAL-ROSETTA-MISSION/ROSINA/RO-C-ROSINA-2-EXT3-V2.0/
             directory DATA/MTP35/COPS/BG (BG = both gauges; NG = nude gauge, RG = ram gauge).
  Local copy /data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2/mtp35_cops_bg/bg_YYYYMMDD_hhmmssmmm_m0326.tab
             (instrument mode M0326): all 1152 files of 2016-09-30 and the last 3 files of 2016-09-29.
  Product    One file = one COPS housekeeping set (object COPS_HK_TABLE, 80-byte ASCII records
             "NAME","FLAG","VALUE","UNIT") for the interval START_TIME..STOP_TIME of the PDS label. The interval is
             60 s at large distance and 4 s in the last part of the descent.
  Fields used
             ROSINA_COPS_PRESSURE_NG, ROSINA_COPS_PRESSURE_RG     gauge pressures [mbar] (N2, 20 C); 0 = no valid value
             ROSINA_COPS_SENSOR_TEMP, ROSINA_COPS_ELECTR_TEMP     sensor and electronics temperature [C] (-> K);
                                                                  the sensor temperature is used as T_RG
             ROSINA_COPS_EMISSION_I_FILAMENT, ..._I_MICROTIP      emission currents [uA] of the NG (filament) and the
                                                                  RG (microtip) ion sources
             ROSINA_COPS_FIL_EMISSION_R_FLAG, ROSINA_COPS_MIC_EMISSION_R_FLAG   emission range (HIGH/LOW)
             ROSINA_COPS_AUTORANGE_NG_FLAG, ROSINA_COPS_AUTORANGE_RG_FLAG       autorange (ON/OFF)
             DATA_QUALITY_ID of the label
  Documentation
             ROSINA Data Users Guide RO-ROS-UG-1000 (DOCUMENT/ROSINA_DATA_USERS_GUIDE.PDF of the data set), sect. 5:
             pressure formulae, n = p/(kT) with T = 293 K, ram gauge flux formula, operation modes.
  Known features of 2016-09-30
             NG emission range HIGH -> LOW at 10:24:40; NG pressure 0 from 10:24:42 to 10:33:14; NG reading constant
             at 2.26e-8 mbar from 10:23:46 to 10:24:38 (saturation of the HIGH range); after 10:33:14 (below ~300 m)
             the NG values alternate between two levels, ~6e-8 and ~9e-8 mbar.
             The srcCG reference table (RosinaMeasurements.h) contains the same pressures with the time of the END
             of each interval.

OUTPUT
  One row per file, sorted by time; the time of a row (utc_mid) is the middle of the interval.
  Columns:
    utc_mid  start  stop  pNG_mbar  pRG_mbar  Tsensor_K  Telectr_K  IemisFIL_uA (NG)  IemisMIC_uA (RG)
    emisRangeFIL  emisRangeMIC  autorangeNG  autorangeRG  quality
  The table data/COPS_HK_20160930.txt is used by LiouvilleCOPS (TemperatureRGMode = hk_table: T_RG of each
  observation interval) and by tools/collect_cops_sources.py (dataset PSA_L2_BG).

usage: cops_hk_table.py [directory with bg_*.tab] [output table]
   e.g. cops_hk_table.py /data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2/mtp35_cops_bg data/COPS_HK_20160930.txt
"""
import glob
import argparse
import os
import sys
from datetime import datetime

from cops_paths import cops_l2_root, prepare_output


def parse(fname):
    hk = {}
    lab = {}
    with open(fname, errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith('"ROSINA_COPS_'):
                fields = [x.strip() for x in line.split('","')]
                name = fields[0].lstrip('"')
                flag = fields[1] if len(fields) > 1 else ""
                value = fields[2] if len(fields) > 2 else ""
                hk[name] = (flag, value)
            elif "=" in line and not line.startswith(" "):
                k, v = line.split("=", 1)
                lab[k.strip()] = v.strip().strip('"')
    return lab, hk


def num(hk, key):
    v = hk.get(key, ("", ""))[1]
    return float(v) if v else float("nan")


def flag(hk, key):
    return hk.get(key, ("", ""))[0] or "NA"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", nargs="?", default=str(cops_l2_root() / "mtp35_cops_bg"),
                        help="raw bg_*.tab directory (default: <COPS_DATA_ROOT>/data_COPS_L2/mtp35_cops_bg)")
    parser.add_argument("output", nargs="?", default="data/COPS_HK_20160930.txt",
                        help="table under cops_output, or absolute path outside LiouvilleCOPS (default: data/COPS_HK_20160930.txt)")
    args = parser.parse_args()

    rows = []
    for fname in sorted(glob.glob(os.path.join(args.directory, "bg_*.tab"))):
        lab, hk = parse(fname)
        if "ROSINA_COPS_SENSOR_TEMP" not in hk:
            sys.exit("Error: %s has no ROSINA_COPS_SENSOR_TEMP" % fname)

        start, stop = lab["START_TIME"], lab["STOP_TIME"]
        rows.append((start, stop,
                     num(hk, "ROSINA_COPS_PRESSURE_NG"), num(hk, "ROSINA_COPS_PRESSURE_RG"),
                     num(hk, "ROSINA_COPS_SENSOR_TEMP") + 273.15, num(hk, "ROSINA_COPS_ELECTR_TEMP") + 273.15,
                     num(hk, "ROSINA_COPS_EMISSION_I_FILAMENT"), num(hk, "ROSINA_COPS_EMISSION_I_MICROTIP"),
                     flag(hk, "ROSINA_COPS_FIL_EMISSION_R_FLAG"), flag(hk, "ROSINA_COPS_MIC_EMISSION_R_FLAG"),
                     flag(hk, "ROSINA_COPS_AUTORANGE_NG_FLAG"), flag(hk, "ROSINA_COPS_AUTORANGE_RG_FLAG"),
                     lab.get("DATA_QUALITY_ID", "NA")))

    rows.sort()
    if not rows:
        sys.exit("Error: no bg_*.tab files found in %s; no output written" % args.directory)
    output = prepare_output(args.output)

    with open(output, "w") as out:
        out.write("# COPS housekeeping from %d L2 bg files in %s\n" % (len(rows), os.path.abspath(args.directory)))
        out.write("# utc_mid start stop pNG_mbar pRG_mbar Tsensor_K Telectr_K IemisFIL_uA IemisMIC_uA emisRangeFIL emisRangeMIC autorangeNG autorangeRG quality\n")
        fmt = "%Y-%m-%dT%H:%M:%S.%f"
        for r in rows:
            t0, t1 = datetime.strptime(r[0], fmt), datetime.strptime(r[1], fmt)
            mid = (t0 + (t1 - t0) / 2).strftime(fmt)[:-3]
            out.write("%s %s %s %.5e %.5e %.3f %.3f %.5e %.5e %s %s %s %s %s\n" % ((mid,) + r))

    print("wrote %s: %d rows, %s .. %s" % (output, len(rows), rows[0][0], rows[-1][1]))


if __name__ == "__main__":
    main()
