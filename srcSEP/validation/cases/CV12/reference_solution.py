#!/usr/bin/env python3
"""Independent closed-system energy identity for CV12.

The reference is deliberately only E_particle+E_wave=constant; it does not
reuse the scatterer's kinematics, so an equal-and-opposite ledger error remains
observable.  The negative control intentionally violates this identity.
"""
import argparse,csv,json,math
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input,encoding="utf-8"));m=1.67262192369e-27;cl=299792458.;v=c["physics"]["particle_speed_m_per_s"];b=v/cl;g=1/math.sqrt(1-b*b);ke=m*cl*cl*g*g*b*b/(g+1)
with open(a.output,"w",newline="",encoding="utf-8") as f:
 w=csv.writer(f);w.writerow(["particle_count","initial_total_energy_j"])
 for n in c["numerics"]["particle_counts"]:w.writerow([n,n*ke+c["physics"]["initial_wave_energy_j"]])
