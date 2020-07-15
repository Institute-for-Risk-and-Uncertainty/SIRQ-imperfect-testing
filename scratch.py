import SIRQModel as Q
import matplotlib.pyplot as plt
import numpy as np

QTest = Q.SIRQ()
PopSize = 6.7e7
Beta = 0.32
gamma = 0.1

DiffSens1e6 = True
PrevSens1e6 = True
PrevSpec1e6 = True
CapSens1e6 = True

if DiffSens1e6:
    QTest.Def_Anti_Test(1, 0, 1e9, 0, 1)
    QTest.Def_Population(0.034*PopSize, 0.01*PopSize, 0.001*PopSize, 0.95*PopSize, 0.004*PopSize, 0.001*PopSize)
    [fig,ax] = plt.subplots(2,2)
    for i, Sens in enumerate([0.5, 0.75, 0.9, 0.98]):
        QTest.Def_Inf_Test(Sens, 0.9, 1e5, 1, 1)
        QTest.Def_Population(0.034*PopSize, 0.01*PopSize, 0.001*PopSize, 0.95*PopSize, 0.004*PopSize, 0.001*PopSize)
        Pop = QTest.Run(0.32, 0.1, Length = 100)
        ax[0,0].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1E6, color = [(3-i)/3,0,i/3])
        ax[1,0].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Susceptible'][0,:], Pop['Recovered'][0,:]), 0)/PopSize, label = r'$\sigma_G={:.2f}$'.format(i), color = [(3-i)/3,0,i/3])
        QTest.Def_Inf_Test(Sens, 0.9, 1.5e5, 1, 1)
        QTest.Def_Population(0.034*PopSize, 0.01*PopSize, 0.001*PopSize, 0.95*PopSize, 0.004*PopSize, 0.001*PopSize)
        Pop = QTest.Run(0.32, 0.1, Length = 100)
        ax[0,1].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, color = [(3-i)/3,0,i/3])
        ax[1,1].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Susceptible'][0,:], Pop['Recovered'][0,:]), 0)/PopSize, label = r'$\sigma_G={:.2f}$'.format(i), color = [(3-i)/3,0,i/3])
    plt.subplots_adjust(hspace = .001, wspace = .001)

if PrevSens1e6:
    QTest.Def_Inf_Test(0.9, 0.9, 0, 0, 1)
    [fig,ax] = plt.subplots(2,4, sharey='row')
    for k, Sens in enumerate([0.5, 0.75, 0.9, 0.98]):
        QTest.Def_Anti_Test(Sens, 0.9, 5e5, 0, 1)
        for p, Prev in enumerate([0.001, 0.003, 0.01, 0.1, 0.25, 0.5]):
            QTest.Def_Population(0.034*PopSize, 0.01*PopSize, 0.001*PopSize, 0.95*PopSize*(1-Prev), 0.004*PopSize, 0.95*PopSize*(Prev))
            Pop = QTest.Run(0.32, 0.1, Length = 365)
            ax[0,k].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, color = [(5-p)/5,0,p/5])
            ax[1,k].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Susceptible'][0,:], Pop['Recovered'][0,:]), 0)/PopSize, label = r'$P_0={:.3f}$'.format(Prev), color = [(5-p)/5,0,p/5])
        ax[1,k].set_xlabel(r'$\sigma_B={:.2f}$'.format(Sens))
    plt.legend()
    plt.subplots_adjust(hspace = .001, wspace = .001)

if PrevSpec1e6:
    QTest.Def_Inf_Test(0.9, 0.9, 0, 0, 1)
    [fig,ax] = plt.subplots(2,4, sharey='row')
    for k, Spec in enumerate([0.5, 0.75, 0.9, 0.98]):
        QTest.Def_Anti_Test(0.9, Spec, 5e5, 0, 1)
        for p, Prev in enumerate([0.001, 0.003, 0.01, 0.1, 0.25, 0.5]):
            QTest.Def_Population(0.034*PopSize, 0.01*PopSize, 0.001*PopSize, 0.95*PopSize*(1-Prev), 0.004*PopSize, 0.95*PopSize*(Prev))
            Pop = QTest.Run(0.32, 0.1, Length = 365)
            ax[0,k].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, color = [(5-p)/5,0,p/5])
            ax[1,k].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Susceptible'][0,:], Pop['Recovered'][0,:]), 0)/PopSize, label = r'$P_0={:.3f}$'.format(Prev), color = [(5-p)/5,0,p/5])
        ax[1,k].set_xlabel(r'$\tau_B={:.2f}$'.format(Spec))
    plt.legend()
    plt.subplots_adjust(hspace = .001, wspace = .001)

if CapSens1e6:
    [fig,ax] = plt.subplots(3,5, sharey='row')
    for i, Cap in enumerate([0, 5e4, 8e4, 1e5, 1.2e5]):
        for k, Rate in enumerate([0.01, 0.05, 0.1]):
            for j, Targ in enumerate([0.05, 0.1, 0.2, 0.3, 0.5, 0.7]):
                QTest.Def_Inf_Test(0.9, 0.9, int(Cap), Targ, 1)
                QTest.Def_Anti_Test(1, 0, int(Rate*PopSize), 0, 7)
                QTest.Def_Population(0.034*PopSize, 0.01*PopSize, 0.001*PopSize, 0.95*PopSize, 0.004*PopSize, 0.001*PopSize)
                Pop = QTest.Run(0.32, 0.1, Length = 365)
                ax[k, i].plot(np.arange(365), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, label = r'$T_0={:.3f}$'.format(Targ), zorder = 100, color = [(5-j)/5,0,j/5])
                for a in [0,1,2]:
                    ax[a, i].plot(np.arange(365), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, color = [0.7,0.7,0.7,1], zorder = 0)
            ax[k, 0].set_ylabel(r'$Rate = {:.2f}$'.format(Rate))
        ax[2,i].set_xlabel(r'$Cap.={:.2e}$'.format(Cap))
    ax[2,0].set_xlabel(r'No Testing')
    plt.subplots_adjust(hspace = .001, wspace = .001)
plt.show()