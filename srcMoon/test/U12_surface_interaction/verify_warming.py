"""Independent quadrature oracle for imposed warming, never observations."""
import csv
import io
import math


def evaluate(stdout,contract):
    from scipy.integrate import quad
    parsed=[row for row in csv.reader(io.StringIO(stdout)) if row and row[0]=='WARMING']
    if len(parsed)!=6:raise ValueError('expected Na/Ar at three warming resolutions')
    comparisons=[];failed=False
    for species in ('Na','Ar'):
        rows=[r for r in parsed if r[1]==species]
        H,start,stop=map(float,rows[0][2:5])
        def hazard(t):
            T=start+(stop-start)*t/H
            if species=='Na':return 1e13*math.exp(-1.85*1.602176634e-19/(1.380649e-23*T))
            return 1.380649e-23*T/6.62607015e-34*math.exp(-28000/(8.31446261815324*T))
        def survival(t):
            integrated,_=quad(hazard,0,t,epsabs=1e-10,epsrel=1e-10)
            return math.exp(-integrated)
        oracle_survival=survival(H)
        oracle_mean,error=quad(survival,0,H,epsabs=1e-9,epsrel=1e-10)
        if error>contract['finest_restricted_residence_error_fraction_H']*H*.01:
            raise RuntimeError('independent quadrature error exceeds its declared budget')
        means=[];survivals=[]
        for row,dt_expected in zip(rows,contract['dt_s']):
            dt,actual_survival,actual_mean=map(float,row[5:8])
            if dt!=dt_expected:raise ValueError('warming timestep differs from frozen contract')
            means.append(actual_mean);survivals.append(actual_survival)
            comparisons.append({'species':species,'dt_s':dt,'model_first_survival':actual_survival,
                'oracle_first_survival':oracle_survival,'model_restricted_residence_s':actual_mean,
                'oracle_restricted_residence_s':oracle_mean,'oracle_error_s':error})
        decreasing=abs(means[2]-means[1])<abs(means[1]-means[0])
        passed=(abs(means[-1]-oracle_mean)<=contract['finest_restricted_residence_error_fraction_H']*H and
                abs(survivals[-1]-oracle_survival)<=contract['finest_survival_absolute_error_limit'] and decreasing)
        failed=failed or not passed
    return {'status':'FAIL' if failed else 'PASS','scope':'U12 production release kernel under imposed warming; not thermal/observational validation','comparisons':comparisons}


def save_artifacts(result,directory):
    """Preserve evaluated analytical rows and show their numerical comparison."""
    directory.mkdir(parents=True,exist_ok=True)
    rows=result['comparisons']
    with (directory/'comparison.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for axis,species in zip(axes,('Na','Ar')):
        group=[r for r in rows if r['species']==species]
        axis.plot([r['dt_s'] for r in group],[r['model_restricted_residence_s'] for r in group],'o-',label='production release kernel')
        axis.axhline(group[0]['oracle_restricted_residence_s'],color='black',label='independent quadrature')
        axis.set(title=species,xlabel='step [s]',ylabel='restricted first residence [s]');axis.legend(fontsize=7)
    fig.suptitle('IMPOSED ANALYTICAL WARMING — NOT OBSERVATIONS')
    for extension in ('png','eps'):fig.savefig(directory/('comparison.'+extension),dpi=150)
    plt.close(fig)
