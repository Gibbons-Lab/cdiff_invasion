import micom
from micom.workflows import workflow
from micom import load_pickle
from micom.media import minimal_medium
from micom.solution import crossover
import pandas as pd
from tqdm import tqdm as progbar
from time import time
from pandas.core.common import SettingWithCopyWarning
import warnings
from sys import argv
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
from os import listdir
import numpy as np

start = time()
logger = micom.logger.logger
try:
    max_procs = snakemake.threads
except NameError:
    max_procs = 10 
 
def calc_rates(args):
    f=args[0]
    frac=args[1]
    smp=f.split('.pickle')[0].split('_')[1]
    cond=f.split('.pickle')[0].split('_')[0]
    try:
        com=load_pickle('%s/%s'%(argv[1],f))
        abund=com.abundances.copy()
        com.taxa = list(abund.index.values)
        sol = com.cooperative_tradeoff(fraction=frac,fluxes=True,pfba=True, rtol=1e-3, atol=1e-4)
        rates = sol.members[["growth_rate",'abundance']].copy()
        rates=rates.drop('medium')
        rates.loc["community",:] = [sol.growth_rate,rates['abundance'].sum()]
        rates['sample'] = smp
        rates['condition']=cond
        rates['tradeoff']= frac
        fluxes = sol.fluxes.copy()
        fluxes['sample'] = smp
        fluxes['condition']=cond
        fluxes['tradeoff']= frac
        return {"rates": rates, "fluxes": fluxes}
    except:
        print('Issue with sample %s'%(smp))
        return None

files = listdir(argv[1])
files = [f for f in files if 'pickle' in f and 'probiotic' in f and 'no-cdiff' in f ]
#files=np.random.choice(files,500,replace=False)
#files = np.random.choice(files,3000)
tradeoffs=[.8]
args=[]
for f in files:
    for t in tradeoffs:
        args.append((f,t))
print('max_procs:%s, args:%s'%(max_procs,len(args)))
res = workflow(calc_rates, args, max_procs)
med=pd.DataFrame()
rates=pd.DataFrame()
fluxes=pd.DataFrame()
i=0
for r in progbar(res):
    if r != None:
        rates = rates.append(r['rates'],sort=False)
        fluxes = fluxes.append(r['fluxes'],sort=False)
    else:
        print('Issue with sample %s'%(args[i][0]))
    i+=1
#Be sure to change 01 to 001 for 1% cdiff samples
rates.to_pickle('{:%m%d%Y}_no-cdiff_growth_rates_probiotic_tradeoff-vanco_pfba.pkl'.format(pd.Timestamp('today')))
fluxes.to_pickle('{:%m%d%Y}_no-cdiff_fluxes_probiotic_tradeoff-vanco_pfba.pkl'.format(pd.Timestamp('today')))
end=time()
print(end-start)
