import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pymc3 as pm3
from tqdm import tqdm

GlobalDeaths = pd.read_csv("time_series_covid19_deaths_global.csv")
UKDeaths = GlobalDeaths.iloc[223].values[4:]

GlobalCases = pd.read_csv('time_series_covid19_confirmed_global.csv')
UKCases = GlobalCases.iloc[223].values[4:]

GlobalRecovered = pd.read_csv('time_series_covid19_recovered_global.csv')
UKRecovered = GlobalRecovered.iloc[223].values[4:]



[figs,ax] = plt.subplots(1,3)
DAx = ax[0]
CAx = ax[1]
RAx = ax[2]

def LogLklhdCOVID(RRate,DRate,InfRate,IncRateSymRate,HosRate,CritRate,Conts,StartIncs,Runs):
Runs = 10000
RecoveryInterval = [0.004,0.004]
DeathInterval = [0.01,0.004]
InfectionInterval = [0.15,0.3]
IncubationInterval = [0.02,0.02]
SymptomaticInterval = [0.6,0.6]
HospitalisedInterval = [0.1,0.1]
CriticalInterval = [0.1,0.1]
InitPop = 6.8e6
contactsInterval = [0,20]
SIRD = np.array([InitPop,0,0,0,0,0,0,0])
SIRD0 = SIRD
SIRD = np.stack((SIRD,[InitPop-12,10,2,0,0,0,0,0]))
SIRD0 = np.stack((SIRD0,))
SIRD = np.array([SIRD]*Runs).astype(int)
SIRD0 = np.array([SIRD0]*Runs).astype(int)
t = 2
Susceptible = SIRD0[:,0,0]
Incubating = SIRD0[:,0,1]
Symptomatic = SIRD0[:,0,2]
Asymptomatic = SIRD0[:,0,3]
Hospitalised = SIRD0[:,0,4]
Critical = SIRD0[:,0,5]
Recoveries = SIRD0[:,0,6]
Deaths = SIRD0[:,0,7]
while True:
    Connections = np.random.binomial(
        np.min((SIRD[:,t-1,0],np.abs(np.random.randint(contactsInterval[0],contactsInterval[1],Runs)*np.sum(SIRD[:,t-1,[2,3]],1))),0),
        SIRD[:,t-1,0]/np.sum(SIRD[:,t-1,0:4],1))

    Incubations = np.random.binomial(Connections,
        np.random.uniform(InfectionInterval[0],InfectionInterval[1],Runs))

    Infections = np.random.binomial(SIRD[:,t-1,1],
        np.random.uniform(IncubationInterval[0],IncubationInterval[1],Runs))

    Symptomatics = np.random.binomial(Infections,
        np.random.uniform(SymptomaticInterval[0],SymptomaticInterval[1],Runs))

    Asymptomatics = Infections - Symptomatics

    Hospitalised = np.random.binomial(SIRD[:,t-1,2],
        np.random.uniform(HospitalisedInterval[0],HospitalisedInterval[1],Runs))

    Criticals = np.random.binomial(SIRD[:,t-1,4],
        np.random.uniform(CriticalInterval[0],CriticalInterval[1],Runs))

    Deaths = np.random.binomial(SIRD[:,t-1,5],
        np.random.uniform(DeathInterval[0],DeathInterval[1],Runs))

    IRecoveries = np.random.binomial(SIRD[:,t-1,1] - Infections,
        np.random.uniform(RecoveryInterval[0],RecoveryInterval[1],Runs))

    SRecoveries = np.random.binomial(SIRD[:,t-1,2] - Hospitalised,
        np.random.uniform(RecoveryInterval[0],RecoveryInterval[1],Runs))

    ARecoveries = np.random.binomial(SIRD[:,t-1,3],
        np.random.uniform(RecoveryInterval[0],RecoveryInterval[1],Runs))

    HRecoveries = np.random.binomial(SIRD[:,t-1,4] - Criticals,
        np.random.uniform(RecoveryInterval[0],RecoveryInterval[1],Runs))

    CRecoveries = np.random.binomial(SIRD[:,t-1,5] - Deaths,
        np.random.uniform(RecoveryInterval[0],RecoveryInterval[1],Runs))

    
    
    SIRD0[:, 0, 1] = SIRD[:, t-1, 1] + Incubations - IRecoveries - Infections
    SIRD0[:, 0, 2] = SIRD[:, t-1, 2] + Symptomatics - SRecoveries - Hospitalised
    SIRD0[:, 0, 3] = SIRD[:, t-1, 3] + Asymptomatics - ARecoveries
    SIRD0[:, 0, 4] = SIRD[:, t-1, 4] + Hospitalised - HRecoveries - Criticals
    SIRD0[:, 0, 5] = SIRD[:, t-1, 5] + Criticals - CRecoveries - Deaths
    SIRD0[:, 0, 6] = SIRD[:, t-1, 6] + (IRecoveries + SRecoveries + ARecoveries
        + HRecoveries + CRecoveries)
    SIRD0[:, 0, 7] = SIRD[:, t-1, 7] + Deaths
    
    SIRD0[:, 0, 0] = InitPop - np.sum(SIRD0[:,0,1:],1)
    SIRD = np.hstack((SIRD,SIRD0))
    t += 1
    if np.all(SIRD[:,t-1,1] == 0) or t == len(UKCases[UKCases>0]):
        break
    #plt.plot(range(0,t),SIRD[:,0],'g')

""" for i in tqdm(range(0,Runs-1,100)):
    CAx.plot(range(0,t),SIRD[i,:,2],':',color=[1,0,1,0.05])
    CAx.plot(range(0,t),SIRD[i,:,1],':',color=[1,0,0,0.1])
    RAx.plot(range(0,t),SIRD[i,:,3],':',color=[0,0,1,0.1])
    DAx.plot(range(0,t),SIRD[i,:,4],':',color=[0,0,0,0.1]) """

CAx.plot(range(0,len(UKCases[UKCases>0])),UKCases[UKCases>0],'k',label='Confirmed Cases')
CAx.plot(range(0,t),np.max(SIRD[:,:,1],0),'b:',label='Incubating')
CAx.plot(range(0,t),np.min(SIRD[:,:,1],0),'b:')
CAx.plot(range(0,t),np.max(SIRD[:,:,2],0),'y:',label='Symptomatic')
CAx.plot(range(0,t),np.min(SIRD[:,:,2],0),'y:')
CAx.plot(range(0,t),np.max(SIRD[:,:,3],0),'c:',label='Asymptomatic')
CAx.plot(range(0,t),np.min(SIRD[:,:,3],0),'c:')
CAx.plot(range(0,t),np.max(SIRD[:,:,4],0),'m:',label='Hospitalised')
CAx.plot(range(0,t),np.min(SIRD[:,:,4],0),'m:')
CAx.plot(range(0,t),np.max(SIRD[:,:,5],0),'r:',label='Critical')
CAx.plot(range(0,t),np.min(SIRD[:,:,5],0),'r:')
CAx.legend()
#plt.plot(range(0,t),SIRD[:,0],'g',label='Susceptible')

RAx.plot(range(0,len(UKCases[UKCases>0])),UKRecovered[UKCases>0],'b',label='Confirmed Recovered')
RAx.plot(range(0,t),np.max(SIRD[:,:,6],0),'b:',label='Recovered')
RAx.plot(range(0,t),np.min(SIRD[:,:,6],0),'b:',label='Recovered')

DAx.plot(range(0,t),np.max(SIRD[:,:,7],0),'k:',label='Deaths')
DAx.plot(range(0,t),np.min(SIRD[:,:,7],0),'k:',label='Deaths')
DAx.plot(range(0,len(UKCases[UKCases>0])),UKDeaths[UKCases>0],'k',label='Confirmed Deaths')

CAx.set_title('Infections')
RAx.set_title('Recovered')
DAx.set_title('Deaths')
DAx.set_yscale('log')
CAx.set_yscale('log')
RAx.set_yscale('log')
#plt.legend()
plt.show()