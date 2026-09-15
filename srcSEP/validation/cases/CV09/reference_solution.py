#!/usr/bin/env python3
"""Independent test-particle planar DSA slope and acceleration-time reference."""
import argparse,csv,json
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input,encoding="utf-8"));x=c["physics"]
with open(a.output,"w",newline="",encoding="utf-8") as f:
 w=csv.writer(f);w.writerow(["compression_ratio","phase_space_index","acceleration_time_s"])
 for r in x["compression_ratios"]:
  u1=x["upstream_speed_m_per_s"];u2=u1/r;t=3/(u1-u2)*(x["kappa_upstream_m2_per_s"]/u1+x["kappa_downstream_m2_per_s"]/u2);w.writerow([r,3*r/(r-1),t])
