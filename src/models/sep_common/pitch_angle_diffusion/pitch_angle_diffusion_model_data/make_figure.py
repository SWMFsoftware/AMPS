from pathlib import Path
import numpy as np
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(13.4,4.3),layout='constrained')
xi=np.linspace(0,5,61)
closure=[]; variance=[]
for x in xi:
 if x==0:
  closure.append(1/3);variance.append(1/3);continue
 Z=2*np.sinh(x)/x;coth=1/np.tanh(x)
 def kernel(mu):
  rho=np.exp(x*mu)/Z
  H=(np.exp(x*mu)*(mu-coth)+np.exp(-x)*(1+coth))/(x*Z)
  return H*H/(.5*(1-mu*mu)*rho)
 variance.append(quad(kernel,-1,1,epsabs=1e-11,epsrel=1e-10)[0])
 closure.append(coth/x-1/(x*x))
axes[0].plot(xi,closure,color='#d97706',label='Diffusion closure, Eq. (112)',lw=2)
axes[0].plot(xi,variance,color='#1261a0',label='Centered variance, Eq. (115)',lw=2)
axes[0].set(xlabel=r'Focusing strength $\xi$',ylabel=r'$\kappa\nu/v^2$',title='Different focusing observables')
axes[0].legend(fontsize=8,frameon=False)
tau=np.linspace(.001,3,301)
exact=(1+np.expm1(-tau)/tau)/3
axes[1].plot(tau,exact,color='#1261a0',label='Integrated pitch-angle trajectory',lw=2)
axes[1].plot(tau,tau/6,color='#d97706',label='Single stationary endpoint',lw=2)
axes[1].set(xlabel=r'Interval $\nu\Delta t$',ylabel=r'$\mathrm{Var}[s]/(2\Delta t)\;\nu/v^2$',title='Angular endpoints do not fix streaming')
axes[1].legend(fontsize=8,frameon=False)
mu=np.linspace(-1,1,301);a=.3;Z=np.log((1+a)/(1-a))/a
axes[2].plot(mu,1/(Z*(1+a*mu)),color='#1261a0',lw=2,label=r'Correct time change: $\rho\propto 1/g$')
axes[2].plot(mu,np.full_like(mu,.5),color='#d97706',lw=2,label=r'Differentiating $gD$: uniform density')
axes[2].set(xlabel=r'Pitch cosine $\mu$',ylabel='Normalized host-time density',title=r'State-dependent clock: $g=1+0.3\mu$')
axes[2].set_xticks([-1,-.5,0,.5,1])
axes[2].legend(fontsize=8,frameon=False,loc='upper right')
for ax in axes:ax.grid(alpha=.2)
fig.suptitle('Reference predictions for the physical implementation contract',fontsize=14)
fig.savefig(root/'physical_contract_checks.png',dpi=180)
