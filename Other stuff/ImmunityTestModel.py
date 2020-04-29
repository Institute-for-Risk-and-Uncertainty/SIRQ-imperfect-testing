import numpy as np
from tqdm import tqdm

def GenCases(Sample,Runs,UKCases):
    # Model Parameters
    [phi, tau, beta, rho, sigma gamma] = Sample

    InitPop = 6.8e6

    #Initial Population Split
    # [Quarantined_S, Susceptible, Infected, Quarantined_I, Recovered]
    InitQuarantined = np.round(InitPop * 0.9)
    InitSusceptible = 0
    InitInfected = np.round((InitPop - InitQuarantined)/2)
    InitRecovered = InitPop - InitQuarantined - InitInfected
    SIRD = np.array([
        InitQuarantined,    # Quarantined_S
        0,                  # Susceptible
        InitInfected,       # Infected
        0,                  # Quarantined_I
        InitRecovered      # Recovered
        ])
    SIRD0 = SIRD
    SIRD = np.array([SIRD]*Runs).astype(int)
    SIRD0 = np.array([SIRD0]*Runs).astype(int)
    t = 2
    while True:
        # How many get tested who are in quarantine?
        QTests = np.min((np.array([TestCap]*Runs),np.random.binomial(SIRD[:, t-1, 0],phi,Runs)),0)
        # Of those, how many aren't actually immune?
        FalseNegs = np.random.binomial(QTests, (1-tau), Runs)
        # How many of those get newly infected by the infected population?
        Infections = np.min((
            SIRD[:, t-1, 1],
            np.random.binomial(SIRD[:, t-1, 2], beta, Runs)
            ), 0)
        # Of the infected population, how many are tested?
        ITests = np.min((np.array([TestCap]*Runs),np.random.binomial(SIRD[:, t-1, 2],rho,Runs)),0)
        # Of those, how many test positive and are quarantined?
        TruePos = np.random.binomial(SIRD[:, t-1, 2], ITests*(1-sigma), Runs)
        # Of the infected population, how many recover?
        IRecoveries = np.random.binomial(SIRD[:, t-1, 2], gamma, Runs)
        # Of the infected, quarantined population how many recover?
        PRecoveries = np.random.binomial(SIRD[:, t-1, 3], gamma, Runs)

        SIRD0[[:, 0, 0]] = SIRD[:, t-1, 0] - FalseNegs
        SIRD0[[:, 0, 1]] = SIRD[:, t-1, 0] + FalseNegs - Infections
        SIRD0[[:, 0, 2]] = SIRD[:, t-1, 0] + Infections - TruePos - IRecoveries
        SIRD0[[:, 0, 3]] = SIRD[:, t-1, 0] + TruePos - PRecoveries
        SIRD0[[:, 0, 4]] = SIRD[:, t-1, 0] + IRecoveries + PRecoveries

        SIRD = np.hstack((SIRD,SIRD0))
        t += 1
        if np.all(SIRD[:,t-1,1] == 0):
            break
    return SIRD, t