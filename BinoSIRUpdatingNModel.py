import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import tikzplotlib

# This is outdated and should not be used for analysis, mostly just here for reference. Shouldn't go in public available version.

seeds = 100
Runs = 100
burnin = 100
Scaling = 0.01
samples = 100
AGTestCap = 10000# (70% of infected, 30% of susceptible)
ABTestCap = 10000
InitPop = 6.7e7
GlobalDeaths = pd.read_csv("time_series_covid19_deaths_global.csv")
UKDeaths = GlobalDeaths.iloc[223].values[4:]

GlobalCases = pd.read_csv('time_series_covid19_confirmed_global.csv')
UKCases = GlobalCases.iloc[223].values[4:]

GlobalRecovered = pd.read_csv('time_series_covid19_recovered_global.csv')
UKRecovered = GlobalRecovered.iloc[223].values[4:]
# [0.2, tau, 0.32, 0.5, sigma, 0.1]
PriorBounds = np.array([
                [0.5, 1],     #phi
                [0.9, 1],     #tauB
                [0.9, 1],    #sigmaB
                [0.5, 1],     #rho 
                [0.9, 1],   #tauG
                [0.9, 1],     #sigmaG
                [0.0, 0.3],   #gamma
                [0.12, 0.32]     #beta
                ])

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

def GenCases(PriorBounds,Runs,ABTestCap,AGTestCap,InitPerc,SISplit=[0,1],TestInt=1, length = 365):
    # Model Parameters
    # [phi, tau, beta, rho, sigma, gamma] = InitParams
   
    InitPop = 6.7e7

    #Initial Population Split
    # [Quarantined_S, Susceptible, Infected, Quarantined_I, Recovered]
    InitQ_S = int(InitPerc[0]*InitPop)
    InitS = int(InitPerc[1]*InitPop)
    InitQ_I = int(InitPerc[2]*InitPop)
    InitI = int(InitPerc[3]*InitPop)
    InitQ_R = int(InitPerc[4]*InitPop)
    InitR = int(InitPerc[5]*InitPop)
    SIRD = np.array([
        InitQ_S,    # Quarantined_S
        InitS,                  # Susceptible
        InitQ_I,       # Infected
        InitI,                  # Quarantined_I
        InitQ_R,      # Recovered
        InitR    # Recovered but Quarantined
        ])
    SIRD0 = SIRD
    SIRD = np.array([np.stack((SIRD,))]*Runs).astype(int)
    SIRD0 = np.array([np.stack((SIRD0,))]*Runs).astype(int)
    t = 1
    #[phi,tauB,sigmaB,rho,tauG,sigmaG,gamma,beta]
    phi = np.random.uniform(PriorBounds[0,0],PriorBounds[0,1],Runs)
    tauB = np.random.uniform(PriorBounds[1,0],PriorBounds[1,1],Runs)
    sigmaB = np.random.uniform(PriorBounds[2,0],PriorBounds[2,1],Runs)
    rho = np.random.uniform(PriorBounds[3,0],PriorBounds[3,1],Runs)
    tauG = np.random.uniform(PriorBounds[4,0],PriorBounds[4,1],Runs)
    sigmaG = np.random.uniform(PriorBounds[5,0],PriorBounds[5,1],Runs)
    gamma = np.random.uniform(PriorBounds[6,0],PriorBounds[6,1],Runs)
    beta = np.random.uniform(PriorBounds[7,0],PriorBounds[7,1],Runs)
    while True:
        

        SIRD = np.hstack((SIRD,SIRD0))

        # How many want to get tested who are in quarantine?
        if t%TestInt==0:
            QTests = np.random.binomial(np.sum((SIRD[:, t-1, [0,-2]]),1),phi,Runs)
        else:
            QTests = np.array([0]*Runs).astype(int)

        # Of the infected population, how many want to get tested?
        ITests = np.random.binomial(np.sum(SIRD[:, t-1, [1, 3]],1), rho, Runs)

        # Apportion tests fairly
        if np.any(ITests > AGTestCap):
            ITests = np.min(([AGTestCap]*Runs,ITests),0).astype(int)
        
        if np.any(QTests > ABTestCap):
            QTests = np.min(([ABTestCap]*Runs,QTests),0).astype(int)

        # Of those, how many aren't actually immune?
        BFalsePos = np.random.binomial(np.nan_to_num(QTests*(SIRD[:, t-1, 0]/np.sum((SIRD[:, t-1, [0, 4]]),1))).astype(int), (1-tauB), Runs)

        BTruePos = np.random.binomial(np.nan_to_num(QTests*(SIRD[:, t-1, 4]/np.sum((SIRD[:, t-1, [0, 4]]),1))).astype(int), sigmaB, Runs)

        # How many of those get newly infected by the infected population?
        # Infections = np.random.binomial(np.min((SIRD[:, t-1, 1],[int(SIRD[j, t-1, 2])*int(SIRD[j, t-1, 1]) for j in range(0,100)]),0).astype(int), beta, Runs)
        Q_I_Beta = 0.0000000*beta
        OpenPop = np.sum(SIRD[:, t-1, [1, 3, 5]],1)+(SIRD[:, t-1, 2]*Q_I_Beta)

        IInfections = np.random.binomial(SIRD[:, t-1, 3], beta*(SIRD[:, t-1, 1]/OpenPop), Runs)

        QIInfections = np.random.binomial(SIRD[:, t-1, 2], Q_I_Beta*(SIRD[:, t-1, 1]/OpenPop), Runs)

        IPopGTests = np.min((SISplit[1]*ITests, SISplit[1]*SIRD[:, t-1, 3]),0).astype(int)

        SPopGTests = np.min((ITests-IPopGTests, SIRD[:, t-1, 1] - IInfections),0).astype(int)
        
        # How many of the tested infected population test positive?
        GTruePos = np.random.binomial(IPopGTests, sigmaG, Runs)

        # How many of the tested susceptible population test positive?
        GFalsePos = np.random.binomial(SPopGTests, (1-tauG), Runs)

        PosPerc = GTruePos/ITests

        # Of the infected population, how many recover?
        IRecoveries = np.random.binomial(SIRD[:, t-1, 3] - GTruePos, gamma, Runs)

        # Of the infected, quarantined population how many recover?
        PRecoveries = np.random.binomial(SIRD[:, t-1, 2], gamma, Runs)

        SIRD0[:, 0, 0] = np.max(([0]*Runs,SIRD[:, t-1, 0] - BFalsePos + GFalsePos),0) # Q_S
        SIRD0[:, 0, 1] = SIRD[:, t-1, 1] + BFalsePos - GFalsePos - IInfections - QIInfections # S
        SIRD0[:, 0, 2] = SIRD[:, t-1, 2] + GTruePos - PRecoveries # Q_I
        SIRD0[:, 0, 3] = SIRD[:, t-1, 3] + IInfections + QIInfections - GTruePos - IRecoveries # I
        SIRD0[:, 0, 4] = SIRD[:, t-1, 4] - BTruePos # Q_R
        SIRD0[:, 0, 5] = SIRD[:, t-1, 5] + IRecoveries + PRecoveries + BTruePos #R
        change = SIRD0-SIRD[:,t-1,:]
        SIRD[:, [t,], :] = SIRD0
        t += 1
        if np.all(SIRD[:,t-1,1] == 0) or t==length:
            break
    return SIRD, t
"""
[fig, ax] = plt.subplots(2,1)
Bounds = np.copy(PriorBounds)
Bounds[1, :] = [0, 0]
Bounds[4, :] = [0.9, 0.9]
InitPerc = [0.95,0.034,0.004,0.01,0.001,0.001]
[SIRD, t] = GenCases(Bounds, Runs, ABTestCap, AGTestCap, InitPerc)
for j in range(0,Runs):
    ax[0].plot(range(0,t), np.sum(SIRD[j, :, [1,3,5]],0), color=[0,0,1,0.1])
    ax[1].plot(range(0, t), np.sum(SIRD[j, :, [2, 3]],0), color =[1,0,0,0.1])
ax[0].set_ylabel('Unquarantined Population')
ax[1].set_ylabel('Infected Population')
plt.subplots_adjust(hspace = .001)


[fig,ax] = plt.subplots(1,5)
Bounds = np.copy(PriorBounds)
Bounds[1, :] = [0.9, 0.9]
Bounds[4, :] = [0.9, 0.9]
for tau in np.linspace(0.5,0.98,5):
    Bounds[1, :] = [tau, tau]
    [SIRD, t] = GenCases(Bounds, Runs, TestCap)
    for i in range(0,5):
        for j in range(0,Runs):
            ax[i].plot(range(0,t),SIRD[j, :, i], color = [1-tau, 0, tau, 0.1])
ax[0].set_xlabel('Q_S')
ax[1].set_xlabel('S')
ax[2].set_xlabel('I')
ax[3].set_xlabel('Q_I')
ax[4].set_xlabel('R')
plt.subplots_adjust(wspace = .001)
fig.suptitle('Red(tau = 0) Blue(tau = 1)')




PriorBounds = np.array([
                [1, 1],     #phi
                [0.9, 0.9],     #tauB
                [0.95, 0.95],    #sigmaB
                [1, 1],     #rho 
                [0.9, 0.9],   #tauG
                [0.9, 0.9],     #sigmaG
                [0.1, 0.1],   #gamma
                [0.32, 0.32]     #beta
                ])
[fig,ax] = plt.subplots(2,2)
Bounds = np.copy(PriorBounds)
Rez = 4
for k, SigmaG in enumerate([0.5, 0.75, 0.9, 0.98]):
    Bounds[5,:] = [SigmaG, SigmaG]
    Bounds[1,:] = [0, 0]
    Bounds[2,:] = [1, 1]
    Prevalence = 0.001
    InitPerc = [0.951*(1-Prevalence),0.034,0.004,0.01,0.951*(Prevalence),0.001]
    [SIRD, t] = GenCases(Bounds, 1, 1e14, 1e5, InitPerc, length = 100)
    ax[1,0].plot(range(0,t),np.sum(SIRD[0, :, [1,3,5]],0)/InitPop)
    ax[0,0].plot(range(0,t),np.sum(SIRD[0, :, [2,3]],0)/1e6, label = r'$\sigma_G={:.2f}$'.format(SigmaG))
    ax[0,0].plot(range(0,t),np.sum(SIRD[0, :, [0,1]],0)/1e6, label = r'$\sigma_G={:.2f}$'.format(SigmaG))
    Prevalence = 0.001
    ax[0,0].plot(range(0,t),np.sum(SIRD[0, :, [4,5]],0)/1e6, label = r'$\sigma_G={:.2f}$'.format(SigmaG))
    Prevalence = 0.001
    Prevalence = 0.001
    InitPerc = [0.951*(1-Prevalence),0.034,0.004,0.01,0.951*(Prevalence),0.001]
    [SIRD, t] = GenCases(Bounds, 1, 1e14, 1.5e5, InitPerc, length = 100)
    ax[1,1].plot(range(0,t),np.sum(SIRD[0, :, [1,3,5]],0)/InitPop)
    ax[0,1].plot(range(0,t),np.sum(SIRD[0, :, [2,3]],0)/1e6)
ax[1,0].set_xlabel('Days - 1e5 TestCap')
ax[1,1].set_xlabel('Days - 1.5e5 TestCap')
plt.subplots_adjust(hspace = .001)
legend = fig.legend()
#fig.suptitle('Red(sigma = 0.5) Blue(sigma = 0.98)')
#tikzplotlib.save('Diffcap1e6.tikz')


[fig,ax] = plt.subplots(2,1)
Bounds = np.copy(PriorBounds)
Bounds[1, :] = [0, 0]
Bounds[4, :] = [0.9, 0.9]
Rez = 5
for k, TestCapi in enumerate(np.linspace(1e5,250000,Rez).astype(int)):
    InitPerc = [0.95,0.034,0.004,0.01,0.001,0.001]
    [SIRD, t] = GenCases(Bounds, Runs, ABTestCap, TestCapi, InitPerc)
    for i in range(0,5):
        for j in range(0,Runs):
            ax[0].plot(range(0,t),np.sum(SIRD[j, :, [1,3,5]],0), color = [(Rez-k)/Rez, 0, k/Rez, 1])
            ax[1].plot(range(0,t),np.sum(SIRD[j, :, [2,3]],0), color = [(Rez-k)/Rez, 0, k/Rez, 1])
    ax[1].plot(range(0,t),np.sum(SIRD[Runs-1, :, [2,3]],0), color = [(Rez-k)/Rez, 0, k/Rez, 1], label = TestCapi)
            
ax[0].set_ylabel('Unquarantined Population')
ax[1].set_ylabel('Infected Population')
plt.subplots_adjust(hspace = .001)
#fig.suptitle('Variable Antigen Test Capacities')
plt.legend()
tikzplotlib.save('AGTestEff.tikz')


[fig,ax] = plt.subplots(2,4,sharey='row')
Rez = 6
for k, Prevalence in enumerate([0.001, 0.003,  0.01, 0.1, 0.25, 0.5]):
    for i, taui in enumerate([0.5,0.75,0.9,0.98]):
        Bounds = np.copy(PriorBounds)
        Bounds[1,:] = [taui,taui]
        InitPerc = [0.951*(1-Prevalence),0.034,0.004,0.01,0.951*(Prevalence),0.001]
        [SIRD, t] = GenCases(Bounds, Runs, ABTestCap*10, AGTestCap, InitPerc)
        ax[1,i].plot(range(0,t),np.sum(SIRD[0, :, [1,3,5]],0)/InitPop)
        CumInf = SIRD[:,1:,:]-SIRD[:,0:-1,:]
        for R in range(0,Runs):
            if i == 3:
                ax[0,i].plot(range(0,t),np.sum(SIRD[R, :, [2, 3]],0)/1e6, label = r'$P_0={:.3f}$'.format(Prevalence))
            else:
                ax[0,i].plot(range(0,t),np.sum(SIRD[R, :, [2, 3]],0)/1e6)
            if k == 0:
                ax[1,i].set_xlabel(r'$\tau={:.2f}$'.format(taui))
plt.subplots_adjust(hspace = .001, wspace = .001)
legend = fig.legend()
tikzplotlib.save('PrevSpec1e6.tikz')

"""
[fig,ax] = plt.subplots(2,4,sharey='row')
Rez = 6
for k, Prevalence in enumerate([0.001, 0.003, 0.01, 0.1, 0.25, 0.5]):
    for i, sigmaABi in enumerate([0.5,0.75,0.9,0.98]):
        Bounds = np.copy(PriorBounds)
        Bounds[2,:] = [sigmaABi, sigmaABi]
        InitPerc = [0.951*(1-Prevalence),0.034,0.004,0.01,0.951*(Prevalence),0.001]
        [SIRD, t] = GenCases(Bounds, Runs, ABTestCap, AGTestCap, InitPerc)
        ax[1,i].plot(range(0,t),np.sum(SIRD[0, :, [1,3,5]],0)/InitPop)
        if i == 3:
            ax[0,i].plot(range(0,t),np.sum(SIRD[0, :, [2,3]],0)/1e6, label = r'$P_0={:.3f}$'.format(Prevalence))
        else:
            ax[0,i].plot(range(0,t),np.sum(SIRD[0, :, [2,3]],0)/1e6)
        if k == 0:
            ax[1,i].set_xlabel(r'$\sigma_B={:.2f}$'.format(sigmaABi))
plt.subplots_adjust(hspace = .001, wspace = .001)
legend = fig.legend()
#tikzplotlib.save('PrevSens1e6.tikz') 
"""
# Certain number of population leaving every day
[figtau,axtau] = plt.subplots(1,4,sharey='row')
plt.subplots_adjust(hspace = .001, wspace = .001)
[figsig,axsig] = plt.subplots(1,4,sharey='row')
plt.subplots_adjust(hspace = .001, wspace = .001)

[figCap,axCap] = plt.subplots(3,5,sharey='row')
plt.subplots_adjust(hspace = .001, wspace = .001)


[figsplit,axsplit] = plt.subplots(3,5,sharey='row', sharex = 'col')
plt.subplots_adjust(hspace = .001, wspace = .001)

Rez = 2
for k, Leaving in enumerate([0.01, 0.05, 0.1]):
    for i, TPrevalence in enumerate([0.05, 0.1, 0.2, 0.5, 0.7, 0.95]):
        for j, Cap in enumerate([0, 5e4, 8e4, 1e5, 1.2e5]):
            Bounds = np.copy(PriorBounds)
            Bounds[1,:] = [0,0]
            Bounds[2,:] = [1,1]
            InitPerc = [0.95,0.034,0.004,0.01,0.001,0.001]
            [SIRD, t] = GenCases(Bounds, Runs, int(Leaving*InitPop), Cap, InitPerc, SISplit = [1-TPrevalence, TPrevalence], TestInt=7)
            for ki in range(0,3):
                if ki != j:
                    axsplit[ki, j].plot(np.arange(0,t,3),np.sum(SIRD[0, 0::3, [2,3]],0)/1e6,
                        color = [
                            0.7,
                            0.7,
                            0.7,
                            1], zorder = 1)
#for k, Leaving in enumerate([0.01, 0.03, 0.04]):
    
    for i, taui in enumerate([0.7,0.8,0.9]):
        Bounds = np.copy(PriorBounds)
        Bounds[1,:] = [0,0]
        Bounds[2,:] = [1,1]
        Bounds[4,:] = [taui, taui]
        InitPerc = [0.95,0.034,0.004,0.01,0.001,0.001]
        [SIRD, t] = GenCases(Bounds, Runs, int(Leaving*InitPop), AGTestCap, InitPerc, TestInt=7)
        if i == 2:
            axtau[i].plot(range(0,t),np.sum(SIRD[0, :, [2,3]],0), label = r'$Rate={:.3f}/day$'.format(Leaving))
        else:
            axtau[i].plot(range(0,t),np.sum(SIRD[0, :, [2,3]],0))
        if k == 0:
            axtau[i].set_xlabel(r'$\tau_B={:.2f}$'.format(taui))
    for i, sigmai in enumerate([0.8,0.9,0.95]):
        Bounds = np.copy(PriorBounds)
        Bounds[1,:] = [0,0]
        Bounds[2,:] = [1,1]
        Bounds[5,:] = [sigmai, sigmai]
        InitPerc = [0.95,0.034,0.004,0.01,0.001,0.001]
        [SIRD, t] = GenCases(Bounds, Runs, int(Leaving*InitPop), AGTestCap*10, InitPerc, TestInt=7, SISplit = [0.3, 0.7])
        if i == 2:
            axsig[i].plot(range(0,t),np.sum(SIRD[0, :, [2,3]],0), label = r'$Rate={:.3f}/day$'.format(Leaving))
        else:
            axsig[i].plot(range(0,t),np.sum(SIRD[0, :, [2,3]],0))
        if k == 0:
            axsig[i].set_xlabel(r'$\sigma_B={:.2f}$'.format(sigmai))
    for i, Cap in enumerate([0, 5e4, 8e4, 1e5, 1.2e5]):
        for l, sigmai in enumerate([0.8,0.9,0.95]):
            for m, split in enumerate([0.3,0.5,0.7]):
                Bounds = np.copy(PriorBounds)
                Bounds[1,:] = [0,0]
                Bounds[2,:] = [1,1]
                #Bounds[4,:] = [taui, taui]
                Bounds[5,:] = [sigmai, sigmai]
                InitPerc = [0.95,0.034,0.004,0.01,0.001,0.001]
                [SIRD, t] = GenCases(Bounds, Runs, int(Leaving*InitPop), Cap, InitPerc, TestInt=7, SISplit=[split,1-split])
                shade = np.min((1,0.5+((sigmai-0.8)/0.15)*0.5))
                if k == 0:
                    axCap[i].set_xlabel(r'$Cap={:.2f}$'.format(Cap))
                if i == 4 and m == 2:
                    axCap[i].plot(range(0,t),np.sum(SIRD[0, :, [2,3]],0),
                        label = r'$Rate={:.2f}/day, \sigma_G={:.2f}$'.format(Leaving, sigmai),
                        color = [
                            ((Rez-k)/Rez)*shade,
                            0,
                            (k/Rez)*shade,
                            1])
                else:
                    axCap[i].plot(range(0,t),np.sum(SIRD[0, :, [2,3]],0),
                        color=[
                            ((Rez-k)/Rez)*shade,
                            0,
                            (k/Rez)*shade,
                            1])
                 
    iRez = 6
    for i, TPrevalence in enumerate([0.05, 0.1, 0.2, 0.3, 0.5, 0.7]):
        for j, Cap in enumerate([0, 5e4, 8e4, 1e5, 1.2e5]):
            Bounds = np.copy(PriorBounds)
            Bounds[1,:] = [0,0]
            Bounds[2,:] = [1,1]

            InitPerc = [0.95,0.034,0.004,0.01,0.001,0.001]
            [SIRD, t] = GenCases(Bounds, Runs, int(Leaving*InitPop), Cap, InitPerc, SISplit = [1-TPrevalence, TPrevalence], TestInt=7)

            CumInf = SIRD[:,1:,:]-SIRD[:,0:-1,:]
            for R in range(0,Runs):
                if j == 4 and k == 2 and R == Runs-1:
                    axsplit[k, j].plot(range(0,t),np.sum(SIRD[R, :, [2, 3]],0)/1e6, label = r'$Prev.={:.3f}$'.format(TPrevalence),
                            color = [
                                ((iRez-i)/iRez),
                                0,
                                (i/iRez),
                                1], zorder = 2)
                else:
                    axsplit[k, j].plot(range(0,t),np.sum(SIRD[R, :, [2, 3]],0)/1e6,
                            color = [
                                ((iRez-i)/iRez),
                                0,
                                (i/iRez),
                                1], zorder = 2)
            for ki in range(0,3):
                if ki != k:
                    axsplit[ki, j].plot(np.arange(0,t,3),np.sum(SIRD[0, 0::3, [2,3]],0)/1e6,
                        color = [
                            0.7,
                            0.7,
                            0.7,
                            1], zorder = 1)
            for Runs in range(1,Runs):
                axsplit[k, j].plot(range(0,t),np.sum(SIRD[Runs, :, [2,3]],0)/1e6,
                        color = [
                            ((iRez-i)/iRez),
                            0,
                            (i/iRez),
                            1], zorder = 2
            
            if k == 2:
                axsplit[k, j].set_xlabel(r'$Cap.={:.2e}$'.format(Cap))
    axsplit[k, 0].set_ylabel(r'$Rate = {:.2f}$'.format(Leaving))
axsplit[2, 0].set_xlabel(r'No Testing')
plt.subplots_adjust(hspace = .001, wspace = .001)

legend = figtau.legend()
legend = figsig.legend()
legend = figCap.legend()

tikzplotlib.save('CapSens1e6.tikz') 
legend = figsplit.legend()
#tikzplotlib.save('PrevSens.tikz') 


[fig,ax] = plt.subplots(1,5)
Bounds = np.copy(PriorBounds)
Bounds[1, :] = [0, 0]
Bounds[4, :] = [0.9, 0.9]
for TestCap in [1e5,2e5,5e5,1e6,5e6]:
    Bounds[4, :] = [sigma, sigma]
    [SIRD, t] = GenCases(Bounds, Runs, TestCap)
    for i in range(0,5):
        for j in range(0,Runs):
            ax[i].plot(range(0,t),SIRD[j, :, i], color = [1-TestCap/5e6, 0, TestCap/5e6, 0.1])
ax[0].set_xlabel('Q_S')
ax[1].set_xlabel('S')
ax[2].set_xlabel('I')
ax[3].set_xlabel('Q_I')
ax[4].set_xlabel('R')
plt.subplots_adjust(wspace = .001)
fig.suptitle('Red(TestCap = 1e5) Blue(TestCap = 5e6)')

gridsplit = 10
Uncert = 1

[fig, ax] = plt.subplots(gridsplit,gridsplit, sharex = 'col', sharey = 'row')
for i, tau in enumerate(np.linspace(PriorBounds[1,0],PriorBounds[1,1],gridsplit)):
    for j, sigma in enumerate(np.linspace(PriorBounds[4,0],PriorBounds[4,1],gridsplit)):
        Bounds = np.copy(PriorBounds)
        #Bounds[1,:] = [tau*0.9,1-((1-tau)*0.9)]
        #Bounds[4,:] = [sigma*0.9,1-((1-sigma)*0.9)]
        Bounds[1, :] = [tau, tau]
        Bounds[4, :] = [sigma, sigma]
        Bounds[:, 0] = Bounds[:, 0] * Uncert
        Bounds[:, 1] = 1 - ((1 - Bounds[:, 1]) * Uncert)
        [SIRD, t] = GenCases(Bounds, Runs, TestCap)
        for k in range(0,Runs):
            ax[i, j].plot(range(0, t), SIRD[k, :, 2], color=[0,0,0,0.1], label='Infected')
            ax[i,j].set_ylim([0,1e7])
        if i == gridsplit - 1:
            ax[i, j].set_xlabel(r'$\sigma={:.2f}$'.format(sigma))
        if j == 0:
            ax[i, j].set_ylabel(r'$\tau={:.2f}$'.format(tau))
plt.subplots_adjust(hspace = .001)
plt.subplots_adjust(wspace = .001)
#plt.tight_layout() """
plt.show()