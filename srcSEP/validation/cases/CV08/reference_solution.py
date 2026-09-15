#!/usr/bin/env python3
"""Independent inverse-Gaussian first-passage CDF for CV08."""
import argparse,csv,json,math
def phi(x):return .5*(1+math.erf(x/math.sqrt(2)))
def cdf(t,L,u,k):
 if t<=0:return 0
 z=math.sqrt(2*k*t);return phi((u*t-L)/z)+math.exp(u*L/k)*phi(-(u*t+L)/z)
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input,encoding="utf-8"));L=c["physics"]["boundary_m"];k=c["physics"]["kappa_m2_per_s"];T=c["numerics"]["maximum_time_s"]
with open(a.output,"w",newline="",encoding="utf-8") as f:
 w=csv.writer(f);w.writerow(["drift_m_per_s","time_s","cdf"])
 for u in c["physics"]["drifts_m_per_s"]:
  for i in range(1,101):w.writerow([u,T*i/100,cdf(T*i/100,L,u,k)])
