#!/usr/bin/env python3
"""Exact shock trajectory and strong-shock DSA reference for IV05."""
import argparse,csv,json
p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--output",required=True);a=p.parse_args();c=json.load(open(a.input));x=c["physics"]
with open(a.output,"w",newline="")as f:
 w=csv.writer(f);w.writerow(["shock_speed_m_per_s","crossing_time_s","q"])
 for v in x["shock_speeds_m_per_s"]:w.writerow([v,(x["node_m"]-x["shock_start_m"])/v,4])
