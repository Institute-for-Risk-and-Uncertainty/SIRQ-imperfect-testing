import SIRQModel as Q
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('font',**{'family':'serif'})
rc('text', usetex=True)

QTest = Q.SIRQ()
PopSize = 6.7e7
Beta = 0.32
gamma = 0.1

DiffSens1e6 = True
PrevSens1e6 = True
PrevSpec1e6 = True
CapSens1e6 = True

if DiffSens1e6:
    QTest.Def_Anti_Test(1, 0, 1e9, 1, 1)
    QTest.Def_Population(0.984*PopSize, 0.01*PopSize, 0.001*PopSize, 0*PopSize, 0.004*PopSize, 0.001*PopSize)
    [fig,ax] = plt.subplots(2,3, sharex = 'col', sharey = 'row', figsize=(5.25102, 5.25102*(4/6)))
    for i, Sens in enumerate([0.5, 0.75, 0.9, 0.98]):
        QTest.Def_Inf_Test(Sens, 0.9, 1e5, 0.8, 1)
        QTest.Def_Population(0.984*PopSize, 0.01*PopSize, 0.001*PopSize, 0*PopSize, 0.004*PopSize, 0.001*PopSize)
        Pop = QTest.Run(0.32, 0.1, Length = 100)
        ax[0,0].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1E6, color = [(3-i)/3,0,i/3])
        ax[1,0].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Susceptible'][0,:], Pop['Recovered'][0,:]), 0)/PopSize, label = r'$\sigma_A={:.2f}$'.format(Sens), color = [(3-i)/3,0,i/3])


        QTest.Def_Inf_Test(Sens, 0.9, 1.5e5, 0.8, 1)
        QTest.Def_Population(0.984*PopSize, 0.01*PopSize, 0.001*PopSize, 0*PopSize, 0.004*PopSize, 0.001*PopSize)
        Pop = QTest.Run(0.32, 0.1, Length = 100)        
        ax[0,1].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, color = [(3-i)/3,0,i/3]) 
        ax[1,1].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Susceptible'][0,:], Pop['Recovered'][0,:]), 0)/PopSize, label = r'$\sigma_A={:.2f}$'.format(Sens), color = [(3-i)/3,0,i/3])

        QTest.Def_Inf_Test(Sens, 0.9, 2e5, 0.8, 1)
        QTest.Def_Population(0.984*PopSize, 0.01*PopSize, 0.001*PopSize, 0*PopSize, 0.004*PopSize, 0.001*PopSize)
        Pop = QTest.Run(0.32, 0.1, Length = 100)
        ax[0,2].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, color = [(3-i)/3,0,i/3])
        ax[1,2].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Susceptible'][0,:], Pop['Recovered'][0,:]), 0)/PopSize, label = r'$\sigma_A={:.2f}$'.format(Sens), color = [(3-i)/3,0,i/3])

    ax[0,0].set_ylabel(r'Infected $\left(10^6\right)$')
    ax[1,0].set_ylabel(r'Exposed Prop.')
    #ax[1,0].set_yticks([0.97,0.98, 0.99,1.0])
    ax[1,0].set_ylim([0.975,1.005])
    ax[1,1].set_xlabel('Days')
    ax[1,1].set_xticks([50,100])
    ax[1,2].set_xticks([50,100])
    ax[0,0].set_xlabel(r'$C_A=1.0\times 10^5$')
    ax[0,1].set_xlabel(r'$C_A=1.5\times 10^5$')
    ax[0,2].set_xlabel(r'$C_A=2.0\times 10^5$')
    ax[0,0].xaxis.set_label_position('top')
    ax[0,1].xaxis.set_label_position('top')
    ax[0,2].xaxis.set_label_position('top')
    plt.subplots_adjust(hspace = .001, wspace = .001)
    plt.legend(bbox_to_anchor = (-0.5,2.15), loc='lower center', ncol = 4)
    plt.gcf().subplots_adjust(bottom=0.15, top = 0.83)
    plt.savefig('Figure5.eps')

if PrevSens1e6:
    QTest.Def_Inf_Test(0.9, 0.9, 0, 0, 1)
    [fig,ax] = plt.subplots(2,4, sharex = 'col', sharey='row', figsize=(5.25102, 5.25102*(4/6)))
    for k, Sens in enumerate([0.5, 0.75, 0.9, 0.98]):
        QTest.Def_Anti_Test(Sens, 0.9, 2e5, 0.8, 1)
        for p, Prev in enumerate([0.001, 0.003, 0.01, 0.1, 0.25, 0.5]):
            QTest.Def_Population(0.034*PopSize, 0.01*PopSize, 0.001*PopSize, 0.95*PopSize*(1-Prev), 0.004*PopSize, 0.95*PopSize*(Prev))
            Pop = QTest.Run(0.32, 0.1, Length = 365)
            ax[0,k].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, color = [(5-p)/5,0,p/5])
            ax[1,k].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Susceptible'][0,:], Pop['Recovered'][0,:]), 0)/PopSize, label = r'$P_0={:.3f}$'.format(Prev), color = [(5-p)/5,0,p/5])
        ax[0,k].set_xlabel(r'$\sigma_B={:.2f}$'.format(Sens))
        ax[0,k].xaxis.set_label_position('top')
    plt.legend()
    ax[1,1].set_xlabel('Days')
    ax[1,1].xaxis.set_label_coords(1, -0.25)
    ax[1,0].set_xticks([0,100,200,300])
    ax[1,0].set_ylim([0, 0.61])
    ax[1,1].set_xticks([0,100,200,300])
    ax[1,2].set_xticks([0,100,200,300])
    ax[1,3].set_xticks([0,100,200,300])
    plt.legend(bbox_to_anchor = (-1,2.15), loc='lower center', ncol = 3)
    plt.gcf().subplots_adjust(bottom=0.15, top = 0.78)
    ax[0,0].set_ylabel(r'Infected $\left(10^6\right)$')
    ax[1,0].set_ylabel(r'Exposed Prop.')
    plt.subplots_adjust(hspace = .001, wspace = .001)
    plt.savefig('Figure6.eps')

if PrevSpec1e6:
    QTest.Def_Inf_Test(0.9, 0.9, 0, 0, 1)
    [fig,ax] = plt.subplots(2,4, sharex = 'col', sharey='row', figsize=(5.25102, 5.25102*(4/6)))
    for k, Spec in enumerate([0.5, 0.75, 0.9, 0.98]):
        QTest.Def_Anti_Test(0.9, Spec, 2e5, 0.8, 1)
        for p, Prev in enumerate([0.001, 0.003, 0.01, 0.1, 0.25, 0.5]):
            QTest.Def_Population(0.034*PopSize, 0.01*PopSize, 0.001*PopSize, 0.95*PopSize*(1-Prev), 0.004*PopSize, 0.95*PopSize*(Prev))
            Pop = QTest.Run(0.32, 0.1, Length = 365)
            ax[0,k].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, color = [(5-p)/5,0,p/5])
            ax[1,k].plot(np.linspace(0, len(Pop['Infectious'][0,:]), len(Pop['Infectious'][0,:])), np.sum((Pop['Infectious'][0,:], Pop['Susceptible'][0,:], Pop['Recovered'][0,:]), 0)/PopSize, label = r'$P_0={:.3f}$'.format(Prev), color = [(5-p)/5,0,p/5])
        ax[0,k].set_xlabel(r'$\tau_B={:.2f}$'.format(Spec))
        ax[0,k].xaxis.set_label_position('top')
    plt.legend()
    ax[1,1].set_xlabel('Days')
    ax[1,1].xaxis.set_label_coords(1, -0.25)
    ax[1,0].set_xticks([0,100,200,300])
    ax[1,1].set_xticks([0,100,200,300])
    ax[1,2].set_xticks([0,100,200,300])
    ax[1,3].set_xticks([0,100,200,300])
    plt.legend(bbox_to_anchor = (-1,2.15), loc='lower center', ncol = 3)
    plt.gcf().subplots_adjust(bottom=0.15, top = 0.78)
    ax[0,0].set_ylabel(r'Infected $\left(10^6\right)$')
    ax[1,0].set_ylabel(r'Exposed Prop.')
    plt.subplots_adjust(hspace = .001, wspace = .001)
    plt.savefig('Figure7.eps')

if CapSens1e6:
    [fig,ax] = plt.subplots(3,5, sharex = 'col', sharey='row', figsize=(5.25102, 5.25102))
    for i, Cap in enumerate([0, 5e4, 1e5, 1.5e5, 2e5]):
        for k, Rate in enumerate([0.01, 0.05, 0.1]):
            for j, Targ in enumerate([0.05, 0.1, 0.2, 0.3, 0.5, 0.7]):
                QTest.Def_Inf_Test(0.9, 0.9, int(Cap), Targ, 1)
                QTest.Def_Anti_Test(1, 0, int(Rate*PopSize), 0, 7)
                QTest.Def_Population(0.034*PopSize, 0.01*PopSize, 0.001*PopSize, 0.95*PopSize, 0.004*PopSize, 0.001*PopSize)
                Pop = QTest.Run(0.32, 0.1, Length = 365)
                ax[k, i].plot(np.arange(365), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, label = r'$T_A={:.3f}$'.format(Targ), zorder = 100, color = [(5-j)/5,0,j/5])
                for a in [0,1,2]:
                    ax[a, i].plot(np.arange(365), np.sum((Pop['Infectious'][0,:], Pop['Q_Infectious'][0,:]), 0)/1e6, color = [0.7,0.7,0.7,1], zorder = 0)
            ax[k, -1].set_ylabel(r'$Rate = {:.2f}$'.format(Rate))
            ax[k, -1].yaxis.set_label_position('right')
        ax[0,i].set_xlabel(r'$C_A={:.1f}$'.format(Cap/1e5))
        ax[0, i].xaxis.set_label_position('top')
    ax[0,0].set_xlabel(r'No Testing')
    ax[1,0].set_ylabel(r'Infected $\left(10^6\right)$')
    ax[2,2].set_xlabel('Days')
    ax[0,2].set_xlabel(r"Test Capacity $\left(10^5\right)$" "\n" r"$C_A$=1.2")
    ax[1,1].xaxis.set_label_coords(1, -0.25)
    plt.legend(bbox_to_anchor = (-1.5,3.27), loc='lower center', ncol = 3)
    plt.gcf().subplots_adjust(top = 0.81)
    plt.subplots_adjust(hspace = .001, wspace = .001)
    plt.savefig('Figure8.eps')
plt.show()