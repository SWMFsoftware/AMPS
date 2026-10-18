#!/usr/bin/env python3
"""Visualize the ROSINA/COPS reference observations embedded in srcCG.

The nude-gauge (NG) and ram-gauge (RG) reference data are compiled into
srcCG/RosinaMeasurements.h as three parallel arrays:
  ObservationTime[nPoints][100], NudeGaugeReferenceData[nPoints], RamGaugeReferenceData[nPoints]
This script parses them directly from the header, so the figure always shows
exactly what the model is fitted against.

Panels:
  1. Full record (2016-09-30 00:00 - 10:39 UTC): NG and RG pressure, log scale
  2. Final descent (4 s cadence part of the record), same quantities
  3. NG/RG ratio with the target ratios hard-coded in RosinaMeasurements_fit.cpp

NG values == 0 are data gaps (not measured zero pressure); they are drawn as
breaks in the line and shaded.

Usage:
  plot_rosina_observations.py [--header PATH] [--out FILE] [--units raw|pa]
                              [--theme light|dark] [--csv FILE] [--show]
"""

import argparse
import os
import re
import sys
from datetime import datetime

import numpy as np
import matplotlib

if "--show" not in sys.argv:
  matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker


#point indices referenced explicitly in RosinaMeasurements_fit.cpp
FIT_POINTS={1005:"1005",1053:"1053",1101:"1101",1149:"1149"}

#NG/RG ratio targets used in RosinaSample::Fit::Fit() (RosinaMeasurements_fit.cpp)
RATIO_TARGETS=[(1101,3.0),(1149,10.0)]

#colors: reference data-viz palette (categorical slots 1-2 + chart chrome)
THEMES={
  "light":{"surface":"#fcfcfb","ink":"#0b0b0b","ink2":"#52514e","muted":"#898781",
           "grid":"#e1e0d9","axis":"#c3c2b7","gap":"#e1e0d9",
           "ng":"#2a78d6","rg":"#eb6834"},
  "dark": {"surface":"#1a1a19","ink":"#f0efec","ink2":"#c3c2b7","muted":"#898781",
           "grid":"#2c2c2a","axis":"#383835","gap":"#2c2c2a",
           "ng":"#3987e5","rg":"#d95926"},
}


def extract_array(text,name):
  """Return the body between the braces of the C array initializer `name`."""
  m=re.search(re.escape(name)+r"\s*(\[[^\]]*\])+\s*=\s*\{",text)
  if m is None:
    sys.exit("Error: array '%s' is not found in the header" % name)

  end=text.find("}",m.end())
  return text[m.end():end]


def load_observations(header):
  with open(header) as f:
    text=f.read()

  #remove comments so that commented-out entries are not parsed
  text=re.sub(r"/\*.*?\*/","",text,flags=re.S)
  text=re.sub(r"//[^\n]*","",text)

  times=re.findall(r'"([^"]+)"',extract_array(text,"ObservationTime"))
  ng=np.array([float(x) for x in re.findall(r"[-+0-9.Ee]+",extract_array(text,"NudeGaugeReferenceData"))])
  rg=np.array([float(x) for x in re.findall(r"[-+0-9.Ee]+",extract_array(text,"RamGaugeReferenceData"))])

  if not (len(times)==len(ng)==len(rg)):
    sys.exit("Error: array lengths differ: time=%i, NG=%i, RG=%i" % (len(times),len(ng),len(rg)))

  t=np.array([datetime.strptime(s,"%Y-%m-%dT%H:%M:%S") for s in times])
  return t,ng,rg


def gap_intervals(mask):
  """Return (first,last) index pairs of consecutive True runs in mask."""
  runs,i,n=[],0,len(mask)

  while i<n:
    if mask[i]:
      j=i
      while j<n and mask[j]: j+=1
      runs.append((i,j-1))
      i=j
    else:
      i+=1

  return runs


def dense_log_ticks(ax,fmt="%.0e"):
  """Label 1-2-5 steps on a log axis that spans less than ~2 decades."""
  ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0,subs=(1.0,2.0,5.0)))
  ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v,pos: fmt % v))
  ax.yaxis.set_minor_locator(mticker.LogLocator(base=10.0,subs=np.arange(1,10)*1.0))
  ax.yaxis.set_minor_formatter(mticker.NullFormatter())


def style_axes(ax,c):
  ax.set_facecolor(c["surface"])
  ax.grid(True,which="major",color=c["grid"],linewidth=0.6)
  ax.set_axisbelow(True)

  for side in ("top","right"):
    ax.spines[side].set_visible(False)
  for side in ("left","bottom"):
    ax.spines[side].set_color(c["axis"])

  ax.tick_params(colors=c["muted"],labelcolor=c["ink2"],labelsize=8)
  ax.yaxis.label.set_color(c["ink2"])
  ax.xaxis.label.set_color(c["ink2"])
  for title in (ax.title,ax._left_title,ax._right_title):
    title.set_color(c["ink"])


def shade_gaps(ax,t,gaps,c,label=True):
  for k,(a,b) in enumerate(gaps):
    #extend the shaded interval to the neighbouring valid samples
    t0=t[max(a-1,0)]
    t1=t[min(b+1,len(t)-1)]
    ax.axvspan(t0,t1,color=c["gap"],alpha=0.7,linewidth=0,zorder=0,
               label="NG data gap" if (label and k==0) else None)


def plot_pressures(ax,t,ng,rg,sel,gaps,c,ylabel,mark_points):
  ng_plot=np.where(ng>0.0,ng,np.nan)

  shade_gaps(ax,t[sel],[g for g in gap_intervals((ng[sel]<=0.0))],c)
  ax.plot(t[sel],ng_plot[sel],color=c["ng"],linewidth=1.3,label="Nude gauge (NG)")
  ax.plot(t[sel],rg[sel],color=c["rg"],linewidth=1.3,label="Ram gauge (RG)")
  ax.set_yscale("log")
  ax.set_ylabel(ylabel)

  if mark_points:
    idx=np.nonzero(sel)[0]
    for i,name in FIT_POINTS.items():
      if i<idx[0] or i>idx[-1]: continue

      #mark on the RG curve (always valid); NG may be a gap at 1005/1053
      ax.plot(t[i],rg[i],marker="o",markersize=5,color=c["rg"],
              markeredgecolor=c["surface"],markeredgewidth=1.2,zorder=5)
      if ng[i]>0.0:
        ax.plot(t[i],ng[i],marker="o",markersize=5,color=c["ng"],
                markeredgecolor=c["surface"],markeredgewidth=1.2,zorder=5)
      ax.annotate(name,(t[i],rg[i]),xytext=(0,-13),textcoords="offset points",
                  ha="center",fontsize=7,color=c["ink2"])

  #direct labels at the right end of each line
  last=np.nonzero(sel)[0][-1]
  for y,txt in ((ng_plot[last],"NG"),(rg[last],"RG")):
    if np.isfinite(y):
      ax.annotate(txt,(t[last],y),xytext=(5,0),textcoords="offset points",
                  va="center",fontsize=8,color=c["ink"],fontweight="bold")


def main():
  here=os.path.dirname(os.path.abspath(__file__))

  p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  p.add_argument("--header",default=os.path.join(here,"..","RosinaMeasurements.h"),
                 help="path to RosinaMeasurements.h (default: ../RosinaMeasurements.h)")
  p.add_argument("--out",default="RosinaObservations.png",help="output image (png, pdf, svg, ...)")
  p.add_argument("--units",choices=["raw","pa"],default="raw",
                 help="raw: values as stored in the header; pa: multiply by 100 (assumes stored values are mbar)")
  p.add_argument("--theme",choices=["light","dark"],default="light")
  p.add_argument("--csv",default=None,help="also write the parsed data as a CSV table")
  p.add_argument("--show",action="store_true",help="open an interactive window")
  args=p.parse_args()

  c=THEMES[args.theme]
  t,ng,rg=load_observations(args.header)
  n=len(t)

  scale,ylabel=(100.0,"Pressure [Pa]") if args.units=="pa" else (1.0,"Pressure [as stored]")
  ng_s,rg_s=ng*scale,rg*scale

  #split the record by cadence: the final descent is where sampling switches to 4 s
  dt=np.array([(t[i+1]-t[i]).total_seconds() for i in range(n-1)])
  fast=np.nonzero(dt<=10.0)[0]
  iDescent=fast[0] if len(fast)>0 else 0

  ratio=np.where((ng>0.0)&(rg>0.0),ng/np.where(rg>0.0,rg,1.0),np.nan)
  gaps=gap_intervals(ng<=0.0)

  #summary to stdout
  print("Header: %s" % os.path.abspath(args.header))
  print("Points: %i,  %s .. %s UTC" % (n,t[0].isoformat(),t[-1].isoformat()))
  print("Cadence: %.0f s for idx 0-%i, %.0f s from idx %i (%s)" % (dt[0],iDescent,dt[-1],iDescent,t[iDescent].time()))
  print("NG gaps (value<=0): %i points in %i intervals" % (int((ng<=0).sum()),len(gaps)))
  for a,b in gaps:
    print("   idx %4i-%4i  %s - %s" % (a,b,t[a].time(),t[b].time()))
  for i,target in RATIO_TARGETS:
    print("NG/RG at point %i: %.2f (fit target > %g)" % (i,ratio[i],target))

  if args.csv is not None:
    with open(args.csv,"w") as f:
      f.write("index,utc,nude_gauge,ram_gauge,ng_over_rg\n")
      for i in range(n):
        f.write("%i,%s,%e,%e,%s\n" % (i,t[i].isoformat(),ng[i],rg[i],"" if np.isnan(ratio[i]) else "%e" % ratio[i]))
    print("CSV written: %s" % args.csv)

  #figure
  plt.rcParams["font.family"]="sans-serif"
  fig,axes=plt.subplots(3,1,figsize=(11,10),gridspec_kw={"height_ratios":[1.1,1.1,0.8]})
  fig.patch.set_facecolor(c["surface"])

  everything=np.ones(n,dtype=bool)
  descent=np.arange(n)>=iDescent

  #panel 1: full record
  ax=axes[0]
  plot_pressures(ax,t,ng_s,rg_s,everything,gaps,c,ylabel,mark_points=False)
  ax.axvline(t[iDescent],color=c["muted"],linewidth=0.8,linestyle="--")
  ax.annotate("cadence 60 s → 4 s",(mdates.date2num(t[iDescent]),1.0),xycoords=("data","axes fraction"),
              xytext=(-4,-4),textcoords="offset points",ha="right",va="top",fontsize=7,color=c["ink2"])
  ax.set_title("ROSINA/COPS reference observations, %s (full record)" % t[0].date(),loc="left",fontsize=11)
  ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))

  #panel 2: final descent
  ax=axes[1]
  plot_pressures(ax,t,ng_s,rg_s,descent,gaps,c,ylabel,mark_points=True)
  ax.set_title("Final descent, %s – %s UTC (points used explicitly by the fit marked)" %
               (t[iDescent].strftime("%H:%M:%S"),t[-1].strftime("%H:%M:%S")),loc="left",fontsize=11)
  ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(0,60,5)))
  dense_log_ticks(ax)

  #panel 3: NG/RG ratio
  ax=axes[2]
  shade_gaps(ax,t,gaps,c,label=False)
  ax.plot(t,ratio,color=c["ink2"],linewidth=1.1)
  ax.set_yscale("log")
  ax.set_ylabel("NG / RG")
  ax.axhline(1.0,color=c["axis"],linewidth=0.8)
  for i,target in RATIO_TARGETS:
    ax.plot(t[i],ratio[i],marker="o",markersize=6,color=c["ink"],markeredgecolor=c["surface"],markeredgewidth=1.2,zorder=5)
    ax.annotate("pt %i: %.1f (fit target >%g)" % (i,ratio[i],target),(t[i],ratio[i]),
                xytext=(-10,0),textcoords="offset points",ha="right",va="center",fontsize=7,color=c["ink2"])
  ax.set_ylim(top=2.0*np.nanmax(ratio))
  dense_log_ticks(ax,fmt="%g")
  ax.set_title("Nude-to-ram gauge ratio",loc="left",fontsize=11)
  ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
  ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))
  ax.set_xlabel("UTC, %s" % t[0].date())

  for ax in axes:
    style_axes(ax,c)

  leg=axes[0].legend(loc="upper left",fontsize=8,frameon=False)
  for txt in leg.get_texts(): txt.set_color(c["ink2"])

  if args.units=="raw":
    fig.text(0.01,0.005,"Values as stored in RosinaMeasurements.h (unit unconfirmed; the Liouville output scales them by 100). "
             "Use --units pa to apply the mbar→Pa factor.",fontsize=7,color=c["muted"])

  fig.tight_layout(rect=(0,0.015,1,1))
  fig.savefig(args.out,dpi=150,facecolor=c["surface"])
  print("Figure written: %s" % args.out)

  if args.show:
    plt.show()


if __name__=="__main__":
  main()
