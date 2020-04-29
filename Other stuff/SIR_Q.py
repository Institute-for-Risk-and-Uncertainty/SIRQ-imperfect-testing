import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#import tikzplotlib

def SIR_Q(sigma, tau, C,days):
    # Total population, N.
    N = 6.7e7

    # days to simulate

    # Initial number of infected and recovered individuals, I0 and R0.
    I0  = 600000/N
    R0  = 1e5/N

    # At start people in quarantine
    QS0  = 0.9
    QI0 = 0

    # Everyone else, S0, is susceptible to infection initially.
    S0 = 1 - I0 - R0 - QS0 - QI0



    # Contact rate (in 1/days).
    beta= 0.32  # for unkown infection people

    # ideal fraction of sick people tested
    pi = 0.8
    # ideal fraction of well people tested
    psi = 1
    # test capacity

    # recovery rate
    gamma = 1/10

    # return from quarentine rate
    chi = 1/14

    # A grid of time points (in days)
    t = np.linspace(0, days, days)

    # The SIR model differential equations.
    def deriv(y, t, N, beta, gamma, sigma, tau, C, pi, psi):
        S, QS, I, QI, R = y

        # calculate fraction to be tested
        elibile_for_test = pi*I + psi*QS

        if elibile_for_test <= C:
            # all that want test get tested
            rho = pi
            phi = psi

        else:
            # some people from both S and I are tested
            rho = pi * C/elibile_for_test
            phi = psi * C/elibile_for_test


        dSdt  = phi*(1-tau)* QS - beta*S*I
        dQSdt = -phi*(1-tau)*QS
        dIdt  = beta*S*I - rho*sigma*I - (1-rho*sigma)*gamma*I
        dQIdt = rho*sigma*I - gamma*QI
        dRdt  = (1-rho*sigma)*gamma*I + gamma*QI

        if 6.7e7*I < 0:
            dIdt = 0

        return dSdt, dQSdt, dIdt, dQIdt, dRdt

    # Initial conditions vector
    y0 = S0, QS0, I0, QI0, R0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, beta, gamma, sigma, tau, C, pi, psi))
    S, QS, I, QI, R = ret.T
    return S, QS, I, QI, R

def plot(S, QS, I, QI, R,days,fname = None):
    t = np.linspace(0, days, days)
    TI = I + QI
    # Plot the data on three separate curves for S(t), I(t) and R(t)
    [fig, ax] = plt.subplots(1,2)
    ax[0].plot(t, S,'b', alpha=0.5, lw=3, label='Susceptible')
    ax[0].plot(t, QS,'c', alpha=0.5, lw=3, label='Quarantined')
    ax[0].set_xlabel('Time(days)')
    ax[0].set_ylabel('Fraction of population Susceptible/Quarantined')


    ax[1].plot(t, (TI),'r', alpha=0.5, lw=3,label='Infected')
    ax[1].set_yscale('log')
    # ax2.set_ylim([1, 10e6])
    ax[1].set_xlabel('Time(days)')
    ax[1].set_ylabel('Fraction of population Infected')


    legend = fig.legend(loc='upper center',ncol=3)

    """ if fname is not None:
        tikzplotlib.save(fname) """
    plt.show()

if __name__ == '__main__':
    days = 100
    S, QS, I, QI, R = SIR_Q(0,0,1,days)
    plot(S, QS, I, QI, R, days, fname = 'SIR_Q_i.tikz')

    days = 200

    S, QS, I, QI, R = SIR_Q(1,0.95,1e5/6.7e7,days)
    plot(S, QS, I, QI, R, days, fname = 'SIR_Q_095.tikz')

    S, QS, I, QI, R = SIR_Q(1,0.65,1e5/6.7e7,days)
    plot(S, QS, I, QI, R, days, fname = 'SIR_Q_065.tikz')

    S, QS, I, QI, R = SIR_Q(1,0.95,1e6/6.7e7,days)
    plot(S, QS, I, QI, R, days, fname = 'SIR_Q_095_1e6.tikz')

    s = 0.8
    for C in [0.5e5,1e5,2.5e5,5e5,1e6]:
        max_dI = []
        ts = np.linspace(0,1,101)
        for t in ts:
            S, QS, I, QI, R = SIR_Q(s,t,C/6.7e7,days)
            max_dI.append(max(6.7e7*max(np.diff(I)),1))
        plt.plot(ts,max_dI,lw = 3,alpha=0.5,label = '%i' %(C/1000))
    plt.xlabel('Specificity')
    plt.ylabel('$\\mathrm{max}(\\dot{I})$')
    plt.yscale('log')
    plt.legend(title = 'Number of daily tests \\\\ (1000s)')
    # tikzplotlib.save('SIR_Q_T.tikz')

    plt.show()
