#!/usr/bin/env python3
"""Independent exact integration of the CV11 exponential rate equation."""
import argparse,csv,json,math
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input,encoding="utf-8"));x=c["physics"];T=c["numerics"]["final_time_s"]
def energy(s,t):
 if s=="growth":z=2*x["growth_rate_per_s"]*t
 elif s=="damping":z=-2*x["damping_rate_per_s"]*t
 elif s=="cancellation":z=0
 else:z=2*(x["sinusoidal_rate_per_s"]/x["angular_frequency_per_s"]*(1-math.cos(x["angular_frequency_per_s"]*t))-x["damping_rate_per_s"]*t)
 return x["initial_energy_j"]*math.exp(z)
with open(a.output,"w",newline="",encoding="utf-8") as f:
 w=csv.writer(f);w.writerow(["scenario","time_s","wave_energy_j"])
 for s in ("growth","damping","cancellation","sign-change"):
  for i in range(101):w.writerow([s,T*i/100,energy(s,T*i/100)])
