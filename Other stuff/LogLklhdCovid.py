import numpy as np
import matplotlib.pyplot as plt

def LogLklhdCOVID([RRate,DRate,InfRate,IncRate,SymRate,HosRate,CritRate,Conts,
    StartIncs],Runs,UKCases,UKRecovered,UKDeaths):
    # Susceptible = SIRD0[:,0,0]
    # Incubating = SIRD0[:,0,1]
    # Symptomatic = SIRD0[:,0,2]
    # Asymptomatic = SIRD0[:,0,3]
    # Hospitalised = SIRD0[:,0,4]
    # Critical = SIRD0[:,0,5]
    # Recoveries = SIRD0[:,0,6]
    # Deaths = SIRD0[:,0,7]

    InitPop = 6.8e6
    SIRD = np.array([InitPop,0,0,0,0,0,0,0])
    SIRD0 = SIRD
    SIRD = np.stack((SIRD,[InitPop-(2+StartIncs),StartIncs,2,0,0,0,0,0]))
    SIRD0 = np.stack((SIRD0,))
    SIRD = np.array([SIRD]*Runs).astype(int)
    SIRD0 = np.array([SIRD0]*Runs).astype(int)
    t = 2
    while True:
        Incubations = np.random.binomial(Conts, InfRate, Runs)

        Infections = np.random.binomial(SIRD[:,t-1,1], IncRate, Runs)

        Symptomatics = np.random.binomial(Infections, SymRate, Runs)

        Asymptomatics = Infections - Symptomatics

        Hospitalised = np.random.binomial(SIRD[:,t-1,2], HosRate, Runs)

        Criticals = np.random.binomial(SIRD[:,t-1,4], CritRate, Runs)

        Deaths = np.random.binomial(SIRD[:,t-1,5], DRate, Runs)

        IRecoveries = np.random.binomial(SIRD[:,t-1,1] - Infections,
            RRate,Runs)

        SRecoveries = np.random.binomial(SIRD[:,t-1,2] - Hospitalised,
            RRate,Runs)

        ARecoveries = np.random.binomial(SIRD[:,t-1,3],
            RRate,Runs)

        HRecoveries = np.random.binomial(SIRD[:,t-1,4] - Criticals,
            RRate,Runs)

        CRecoveries = np.random.binomial(SIRD[:,t-1,5] - Deaths,
            RRate,Runs)

        
        
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
    
    LogLklhd = 0
    for tt = range(1,t):
        CLogLklhd = np.log(np.sum(([SIRD[ii,tt,2:6] - UKCases[tt] for ii in range(0,Runs)]*(1/Runs))))
        DLogLklhd = np.log(np.sum(([SIRD[ii,tt,2:6] - UKCases[tt] for ii in range(0,Runs)]*(1/Runs))))
        RLogLklhd = np.log(np.sum(([SIRD[ii,tt,2:6] - UKCases[tt] for ii in range(0,Runs)]*(1/Runs))))
        LogLklhd += (CLogLklhd + DLogLklhd + RLogLklhd)

    return LogLklhd