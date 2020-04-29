import numpy as np
import matplotlib.pyplot as plt

def LogLklhdCOVIDSimple([RRate,DRate,InfRate,IncRate,SymRate,HosRate,CritRate,
    Conts,StartIncs],Runs,UKCases,UKRecovered,UKDeaths):
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
        Cases = np.random.binomial(Conts, InfRate, Runs)
        Recoveries = np.random.binomial(Conts, InfRate, Runs)
        Deaths = np.random.binomial(Conts, InfRate, Runs)

        SIRD0[:, 0, 1] = SIRD[:, t-1, 1] + Cases - Recoveries - Deaths
        SIRD0[:, 0, 2] = SIRD[:, t-1, 2] + Recoveries
        SIRD0[:, 0, 3] = SIRD[:, t-1, 3] + Deaths
        SIRD0[:, 0, 0] = InitPop - np.sum(SIRD0[:,0,1:],1)
        SIRD = np.hstack((SIRD,SIRD0))
        t += 1
        if np.all(SIRD[:,t-1,1] == 0) or t == len(UKCases[UKCases>0]):
            break
    
    LogLklhd = 0
    for tt = range(1,t):
        CLogLklhd = np.log(np.sum(([SIRD[ii,tt,1] - UKCases[tt] for ii in range(0,Runs)]*(1/Runs))))
        DLogLklhd = np.log(np.sum(([SIRD[ii,tt,2] - UKCRecovered[tt] for ii in range(0,Runs)]*(1/Runs))))
        RLogLklhd = np.log(np.sum(([SIRD[ii,tt,3] - UKDeaths[tt] for ii in range(0,Runs)]*(1/Runs))))
        LogLklhd += (CLogLklhd + DLogLklhd + RLogLklhd)

    return LogLklhd