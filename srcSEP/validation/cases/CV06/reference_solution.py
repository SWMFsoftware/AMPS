#!/usr/bin/env python3
"""Independent eigenfunction reference for CV06.

For D_mumu=D0(1-mu^2), the pitch-angle operator has Legendre eigenfunctions:
``a_l(t)=a_l(0) exp[-l(l+1)D0 t]``.  This script does not import model code.
"""
import argparse, csv, json, math
p=argparse.ArgumentParser(); p.add_argument("--input",required=True); p.add_argument("--output",required=True); a=p.parse_args()
c=json.load(open(a.input,encoding="utf-8")); d0=c["physics"]["d0_per_s"]; eps=c["physics"]["initial_perturbation"]
with open(a.output,"w",newline="",encoding="utf-8") as f:
 w=csv.writer(f); w.writerow(["mode","time_s","coefficient"])
 for ell in c["physics"]["modes"]:
  for i in range(c["numerics"]["sample_count"]+1):
   t=c["numerics"]["final_time_s"]*i/c["numerics"]["sample_count"]
   w.writerow([ell,t,eps/(2*ell+1)*math.exp(-ell*(ell+1)*d0*t)])
