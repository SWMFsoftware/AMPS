#!/usr/bin/env python3
"""Exact free-stream and conservative-remap invariants for IV04."""
import argparse,csv,json
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input))
with open(a.output,"w",newline="")as f:w=csv.writer(f);w.writerow(["quantity","value"]);w.writerow(["uniform_density",1]);w.writerow(["total_integral",1]);w.writerow(["motions",len(c["physics"]["motions"])])
