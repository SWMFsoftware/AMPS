#!/usr/bin/env python3
"""Linear early-time growth and closed-energy reference for IV06."""
import argparse,csv,json,math
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input));x=c["physics"];dt=c["numerics"]["dt_s"]
with open(a.output,"w",newline="")as f:w=csv.writer(f);w.writerow(["quantity","value"]);w.writerow(["first_step_wave_j",x["initial_wave_j"]*math.exp(2*x["growth_per_s"]*.9*dt)]);w.writerow(["resonant_bin",3]);w.writerow(["total_energy_j",101])
