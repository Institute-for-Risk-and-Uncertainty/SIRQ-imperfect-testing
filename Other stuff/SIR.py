import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
​# Total population, N.
N = 10000
​
# days to simulate
days = 200
​
# Initial number of infected and recovered individuals, I0 and R0.
I0  = 50
R0  = 0
​
# At start noone has been tested
Tp0 = 0
Tn0 = 0
​
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0
​
# Contact rate (in 1/days).
## How many infections people do
beta1 = 0.8   # for unkown infection people
beta2 = 0.01  # for people that tested positive
beta3 = 0.5   # for people that tested regative ##Assume same as unkown
​
# fraction of infected people tested
rho = 0.2
​
# sensitivity of test
sigma = 0.99
​
# recovery rate
gamma1 = 1/10
gamma2 = 1/10
gamma3 = 1/10
​
# A grid of time points (in days)
t = np.linspace(0, days, days)
​
# The SIR model differential equations.
def deriv(y, t, N, beta1, beta2, beta3, rho, gamma1, gamma2, gamma3):
    S, I, Tp, Tn, R = y
    dSdt  = -beta1 * S * I / N - beta2 * S * Tp / N - beta3 * S * Tn / N
    dIdt  = beta1 * S * I / N + beta2 * S * Tp / N + beta3 * S * Tn / N - rho*I - gamma1*(1-rho)*I
    dTpdt = sigma*rho*I - gamma2*Tp
    dTndt = (1-sigma)*rho*I - gamma3*Tn
    dRdt  = gamma1*(1-rho)*I + gamma2*Tp + gamma3*Tn
    return dSdt, dIdt, dTpdt, dTndt, dRdt
​
# Initial conditions vector
y0 = S0, I0, Tp0, Tn0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta1, beta2, beta3, rho, gamma1, gamma2, gamma3))
S, I, Tp, Tn, R = ret.T
print(R[-1]/N)
# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t, S/N,'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, (I+Tp+Tn)/N,'r', alpha=0.5, lw=2, label='Total Infections')
ax.plot(t, I/N,'m', alpha=0.5, lw=2, label='Unkown Infections')
ax.plot(t, Tp/N,'y', alpha=0.5, lw=2, label='Tested Positive')
ax.plot(t, Tn/N,'c', alpha=0.5, lw=2, label='Tested Negative')
ax.plot(t, R/N,'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
ax.set_ylim(0,1)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
# for spine in ('top', 'right', 'bottom', 'left'):
#     ax.spines[spine].set_visible(False)
plt.show()