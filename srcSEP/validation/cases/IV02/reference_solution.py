#!/usr/bin/env python3
"""Independent zero-flux stationary pitch-angle distribution for IV02."""
import argparse,csv,json,math
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input));bins=40
with open(a.output,"w",newline="")as f:
 w=csv.writer(f);w.writerow(["ratio","bin","mu_left","mu_right","probability","mean_mu"])
 for q in c["physics"]["focusing_ratios"]:
  z=2*math.sinh(q)/q if q else 2.;mean=1/math.tanh(q)-1/q if q else 0
  for i in range(bins):
   l=-1+2*i/bins;r=-1+2*(i+1)/bins;prob=(math.exp(q*r)-math.exp(q*l))/(q*z) if q else 1/bins;w.writerow([q,i,l,r,prob,mean])
