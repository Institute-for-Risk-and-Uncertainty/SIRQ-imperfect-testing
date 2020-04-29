import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pymc3 as pm3
from tqdm import tqdm

seeds = 10
Runs = 40
burnin = 1000
Scaling = 0.00001
samples = 100

GlobalDeaths = pd.read_csv("time_series_covid19_deaths_global.csv")
UKDeaths = GlobalDeaths.iloc[223].values[4:]

GlobalCases = pd.read_csv('time_series_covid19_confirmed_global.csv')
UKCases = GlobalCases.iloc[223].values[4:]

GlobalRecovered = pd.read_csv('time_series_covid19_recovered_global.csv')
UKRecovered = GlobalRecovered.iloc[223].values[4:]

#[RRate,DRate,InfRate,IncRate,SymRate,HosRate,CritRate,Conts,StartIncs]
PriorBounds = np.array([
                [0, 0.1],     #RRate
                [0, 0.1],     #DRate
                [0, 0.1],     #InfRate
                [0, 1],     #IncRate
                [0, 1],     #SymRate
                [0, 1],     #HosRate
                [0, 1],     #CritRate
                [0, 20],    #Conts
                [0, 1000]   #StartIncs
                ])

VarsinUse = [0,1,2,7]
Seed = np.random.uniform(
    low = PriorBounds[VarsinUse, 0],
    high = PriorBounds[VarsinUse, 1],
    size = (seeds, 4))
Seed[:,3] = np.random.randint(
    PriorBounds[7, 0],
    PriorBounds[7, 1],
    size = (seeds))
""" Seed[:,8] = np.random.randint(
    PriorBounds[8, 0],
    PriorBounds[8, 1],
    size = (100)) """

if seeds>1:
    PropCov = np.cov(Seed.T)
else:
    PropCov = np.eye(4)*Seed**0.5

PriorDens = np.sum(np.log(1/(np.ptp(PriorBounds[VarsinUse,:],1))))
def PropGen(Seed,PropCov,Bounds):
    Bounded = np.array([False]*len(Seed))
    Prop = np.copy(Seed)
    while True:
        Prop[Bounded == False] = [
            np.random.multivariate_normal(Seedi,PropCov) for 
            Seedi in Seed[Bounded == False]]
        Prop[Bounded == False,-1] = np.round(Prop[Bounded == False,-1])
        Bounded[Bounded == False] = [
            np.all(Propi>=Bounds[:,0]) and 
            np.all(Propi<=Bounds[:,1]) for 
            Propi in Prop[Bounded == False,:]]
        if np.all(Bounded):
            break
    return Prop

def GenCases(Sample,Runs,UKCases):
    [InfRate, RecRate, DRate, Conts] = Sample
    # Susceptible = SIRD0[:,0,0]
    # Cases = SIRD0[:,0,1]
    # Recovered = SIRD0[:,0,2]
    # Deaths = SIRD0[:,0,3]
    InitPop = 6.8e6
    SIRD = np.array([InitPop,0,0,0])
    SIRD0 = SIRD
    SIRD = np.stack((SIRD,[InitPop-2,2,0,0]))
    SIRD0 = np.stack((SIRD0,))
    SIRD = np.array([SIRD]*Runs).astype(int)
    SIRD0 = np.array([SIRD0]*Runs).astype(int)
    t = 2
    while True:
        Cases = np.random.binomial(np.min((SIRD[:,t-1,0],SIRD[:,t-1,1]*Conts),0).astype(int), InfRate, Runs)
        Recoveries = np.random.binomial(SIRD[:, t-1, 1], RecRate, Runs)
        Deaths = np.random.binomial(SIRD[:, t-1, 1] - Recoveries, DRate, Runs)
        SIRD0[:, 0, 1] = SIRD[:, t-1, 1] + Cases - Recoveries - Deaths
        SIRD0[:, 0, 2] = SIRD[:, t-1, 2] + Recoveries
        SIRD0[:, 0, 3] = SIRD[:, t-1, 3] + Deaths
        SIRD0[:, 0, 0] = InitPop - np.sum(SIRD0[:,0,1:],1)
        SIRD = np.hstack((SIRD,SIRD0))
        t += 1
        if np.all(SIRD[:,t-1,1] == 0) or t == len(UKCases[UKCases>0]):
            break
    return SIRD, t

def LogLklhdCOVIDSimple(Sample,Runs,UKCases,UKRecovered,UKDeaths):
    [SIRD, t] = GenCases(Sample,Runs,UKCases)
    LogLklhd = 0
    for tt in range(1,t):
        CLogLklhd = np.nan_to_num(np.log(np.sum(([(np.abs(SIRD[ii,tt,1] - UKCases[tt])*(1/Runs)) for ii in range(0,Runs)]))),neginf=0, posinf=0)
        DLogLklhd = np.nan_to_num(np.log(np.sum(([(np.abs(SIRD[ii,tt,2] - UKRecovered[tt])*(1/Runs)) for ii in range(0,Runs)]))),neginf=0, posinf=0)
        RLogLklhd = np.nan_to_num(np.log(np.sum(([(np.abs(SIRD[ii,tt,3] - UKDeaths[tt])*(1/Runs)) for ii in range(0,Runs)]))),neginf=0, posinf=0)
        LogLklhd += (CLogLklhd + DLogLklhd + RLogLklhd)

    return LogLklhd

MHSamples = [np.zeros(np.shape(Seed))]*samples
print('Calculating Seed Likelihood...')
SeedLklhd = np.array([LogLklhdCOVIDSimple(Sample,Runs,UKCases,UKRecovered,UKDeaths) for Sample in tqdm(Seed)])
AccRatio = 0
print('Metropolis Hastings Samples being generated...')
for i in tqdm(range(0,burnin + samples)):
    Prop = PropGen(Seed,PropCov*Scaling,PriorBounds[VarsinUse])
    PropLklhd = np.array([LogLklhdCOVIDSimple(Sample,Runs,UKCases,UKRecovered,UKDeaths) for Sample in Prop])
    AccSamples = np.array([(SeedLklhd[ii] - PropLklhd[ii])>=np.log(np.random.rand()) for ii, Propi in enumerate(PropLklhd)])
    Seed[AccSamples == True] = Prop[AccSamples == True]
    SeedLklhd[AccSamples == True] = np.copy(PropLklhd[AccSamples == True])
    if i >= burnin:
        MHSamples[i-burnin] = np.copy(Seed)
        AccRatio += sum(AccSamples)
    PropCov = np.cov(Seed.T)
MHSamples = np.vstack(MHSamples)
print(AccRatio/(np.shape(MHSamples)[0]))
[fig,ax] = plt.subplots(len(VarsinUse),len(VarsinUse),)
for i in range(0,len(VarsinUse)):
    for j in range(0,len(VarsinUse)):
        if i != j:
            ax[i,j].scatter(MHSamples[:,i],MHSamples[:,j],s=2)
        else:
            ax[i,j].hist(MHSamples[:,i])
ax[3,0].set_xlabel = 'Infection Rate'
ax[3,1].set_xlabel = 'Recovery Rate'
ax[3,2].set_xlabel = 'Death Rate'
ax[3,3].set_xlabel = 'Connections'
ax[0,3].set_xlabel = 'Infection Rate'
ax[1,3].set_xlabel = 'Recovery Rate'
ax[2,3].set_xlabel = 'Death Rate'
ax[3,3].set_xlabel = 'Connections'

CMax = np.array(np.zeros(len(UKCases[UKCases>0])))
CMin = np.array(np.zeros(len(UKCases[UKCases>0])))
RMax = np.array(np.zeros(len(UKCases[UKCases>0])))
RMin = np.array(np.zeros(len(UKCases[UKCases>0])))
DMax = np.array(np.zeros(len(UKCases[UKCases>0])))
DMin = np.array(np.zeros(len(UKCases[UKCases>0])))

print('Calculating Bounds...')
for Sample in tqdm(MHSamples):
    [SIRD, t] = GenCases(Sample,Runs,UKCases)
    CMax[0:t] = [np.max((CMax[tt],np.max(SIRD[:,tt,1]))) for tt in range(0,t)]
    RMax[0:t] = [np.max((CMax[tt],np.max(SIRD[:,tt,2]))) for tt in range(0,t)]
    DMax[0:t] = [np.max((CMax[tt],np.max(SIRD[:,tt,3]))) for tt in range(0,t)]
    CMin[0:t] = [np.min((CMin[tt],np.min(SIRD[:,tt,1]))) for tt in range(0,t)]
    RMin[0:t] = [np.min((CMin[tt],np.min(SIRD[:,tt,2]))) for tt in range(0,t)]
    DMin[0:t] = [np.min((CMin[tt],np.min(SIRD[:,tt,3]))) for tt in range(0,t)]

[fig,ax] = plt.subplots(1,3)
ax[0].plot(range(0,len(UKCases[UKCases>0])),UKCases[UKCases>0],'m',range(0,len(UKCases[UKCases>0])),CMax,'m:',range(0,len(UKCases[UKCases>0])),CMin,'m:')
ax[1].plot(range(0,len(UKCases[UKCases>0])),UKRecovered[UKCases>0],'b',range(0,len(UKCases[UKCases>0])),RMax,'b:',range(0,len(UKCases[UKCases>0])),RMin,'b:')
ax[2].plot(range(0,len(UKCases[UKCases>0])),UKDeaths[UKCases>0],'k',range(0,len(UKCases[UKCases>0])),DMax,'k:',range(0,len(UKCases[UKCases>0])),DMin,'k:')
plt.show()
np.savetxt('MCSamples.csv', MHSamples, delimiter=',')