#!/usr/bin/env python3
"""Independent cell-integrated sinusoidal characteristic for CV10."""
import argparse,csv,json,math
def avg(i,n,s):
 dx=1/n;l=i*dx-s;r=(i+1)*dx-s;return 1+.25*(math.cos(2*math.pi*l)-math.cos(2*math.pi*r))/(2*math.pi*dx)
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input,encoding="utf-8"));n=max(c["numerics"]["resolutions"]);t=c["numerics"]["duration_s"]
with open(a.output,"w",newline="",encoding="utf-8") as f:
 w=csv.writer(f);w.writerow(["branch","cell","wave_action_j"])
 for branch,sign in (("plus",1),("minus",-1)):
  for i in range(n):w.writerow([branch,i,avg(i,n,sign*t)/n])
