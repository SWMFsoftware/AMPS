#!/usr/bin/env python3
"""Independent high-resolution Parker-spiral characteristic for IV01."""
import argparse,csv,json,math
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input));x=c["physics"];omega=2.86533e-6;r0=x["inner_radius_m"];r1=x["outer_radius_m"];v=x["particle_speed_m_per_s"];mu0=.65
with open(a.output,"w",newline="")as f:
 w=csv.writer(f);w.writerow(["wind_m_per_s","path_length_m","arrival_time_s","final_mu"])
 for wind in x["wind_speeds_m_per_s"]:
  aa=omega/wind
  def B(r):return (r0/r)**2*math.sqrt(1+(aa*r)**2)/math.sqrt(1+(aa*r0)**2)
  n=200000;dr=(r1-r0)/n;length=time=0
  for i in range(n):
   r=r0+(i+.5)*dr;ds=math.sqrt(1+(aa*r)**2)*dr;mu=math.sqrt(max(0,1-(1-mu0*mu0)*B(r)));length+=ds;time+=ds/(v*mu)
  w.writerow([wind,length,time,math.sqrt(max(0,1-(1-mu0*mu0)*B(r1)))])
