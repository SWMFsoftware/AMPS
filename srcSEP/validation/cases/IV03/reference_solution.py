#!/usr/bin/env python3
"""Manifest the independent analytic manufactured-solution definition."""
import argparse,csv,json
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input))
with open(a.output,"w",newline="")as f:w=csv.writer(f);w.writerow(["quantity","value"]);w.writerow(["design_order",2]);w.writerow(["positive",1]);w.writerow(["operators",len(c["physics"]["operators"])])
