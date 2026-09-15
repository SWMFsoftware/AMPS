#!/usr/bin/env python3
"""Independent moment/front reference for the symmetric telegraph process."""
import argparse,csv,json,math
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input,encoding="utf-8"));v=c["physics"]["speed_m_per_s"]
with open(a.output,"w",newline="",encoding="utf-8") as f:
 w=csv.writer(f);w.writerow(["rate_per_s","normalized_time","front_mass","msd_m2","effective_kappa_m2_per_s"])
 for nu in c["physics"]["switching_rates_per_s"]:
  for q in c["numerics"]["normalized_times"]:
   t=q/nu;msd=v*v*(t/nu-(1-math.exp(-2*nu*t))/(2*nu*nu));w.writerow([nu,q,math.exp(-q),msd,msd/(2*t)])
