import numpy as np

class SIRQ():
    def __init__(self, Inf_Test = None, Anti_Test = None, Population = None):
        self.Inf_Test = Inf_Test
        self.Anti_Test = Anti_Test
        self.Population = Population
    
    def Def_Anti_Test(self, Sensitivity, Specificity, Capacity, Targeting = 0, Interval = 1):
        self.Anti_Test = {
            'Sensitivity': Sensitivity,
            'Specificity': Specificity,
            'Capacity': Capacity,
            'Targeting': Targeting,
            'Interval': Interval
            }
    
    def Def_Inf_Test(self, Sensitivity, Specificity, Capacity,  Targeting = 0, Interval = 1):
        self.Inf_Test = {
            'Sensitivity': Sensitivity,
            'Specificity': Specificity,
            'Capacity': Capacity,
            'Targeting': Targeting,
            'Interval': Interval
            }
    
    def Def_Population(self, Susceptible, Infectious, Recovered, Q_Susceptible, Q_Infectious, Q_Recovered):
        self.Population = {
            'Susceptible': [int(Susceptible),],
            'Infectious': [int(Infectious),],
            'Recovered': [int(Recovered),],
            'Q_Susceptible': [int(Q_Susceptible),],
            'Q_Infectious': [int(Q_Infectious),],
            'Q_Recovered': [int(Q_Recovered),]
            }

    def Run(self, Beta, gamma, Runs = 1, Inf_Test = None, Anti_Test = None, Population = None, Length = None):
        if Inf_Test is None:
            Inf_Test = self.Inf_Test
        if Anti_Test is None:
            Anti_Test = self.Anti_Test
        if Population is None:
            Population = self.Population

        assert all((not Inf_Test is None, not Anti_Test is None, not Population is None)), 'Finish Defining Problem'
        
        for k in Population.keys():
            Population[k] = np.array([Population[k]*Runs]).T

        t = 2
        while True:
            S_Quarantined = 0
            I_Quarantined = 0
            if t%Inf_Test['Interval'] == 0:
                if Inf_Test['Targeting'] > 0:
                    Mistargeted_Exposed = np.min((Population['Susceptible'][:,-1], np.max(([0]*Runs,Population['Infectious'][:,-1]/Inf_Test['Targeting']-Population['Infectious'][:,-1]), 0)), 0).astype(int)
                else:
                    Mistargeted_Exposed =  Population['Susceptible'][:,-1]

                Inf_Cap = np.min(np.vstack(([Inf_Test['Capacity']]*Runs, Mistargeted_Exposed + Population['Infectious'][:,-1])),0).astype(int)

                S_Tested = np.array([np.random.hypergeometric(Mistargeted_Exposed[i], Population['Infectious'][i,-1], I) if I >0 else 0 for i, I in enumerate(Inf_Cap)])
                S_Quarantined = np.random.binomial(S_Tested, 1 - Inf_Test['Specificity'])

                I_Tested = Inf_Cap - S_Tested
                I_Quarantined = np.random.binomial(I_Tested, Inf_Test['Sensitivity'])
            Exposed_Population = np.array(np.sum([Population[K][:,-1] for K in ['Susceptible', 'Infectious', 'Recovered']],0), dtype=int) - S_Quarantined - I_Quarantined
            Exposed_Prevalence = (Population['Infectious'][:,-1] - I_Quarantined)/Exposed_Population
            S_Infected = np.min(np.vstack((Population['Susceptible'][:, -1] - S_Quarantined,np.random.binomial(Population['Infectious'][:,-1], Beta*((Population['Susceptible'][:, -1] - S_Quarantined)/Exposed_Population)))), 0).astype(int)
            
            
            if t%Anti_Test['Interval'] == 0 and Anti_Test['Capacity']>0:
                if Anti_Test['Targeting'] > 0:
                    Mistargeted_Quarantined = np.min((Population['Q_Susceptible'][:,-1], np.max(([0]*Runs,Population['Q_Recovered'][:,-1]/Anti_Test['Targeting']-Population['Q_Recovered'][:,-1]), 0)), 0).astype(int)
                else:
                    Mistargeted_Quarantined = Population['Q_Susceptible'][:,-1]

                Anti_Cap = np.min(np.vstack(([Anti_Test['Capacity']]*Runs, Mistargeted_Quarantined + Population['Q_Recovered'][:,-1])),0).astype(int)
            
                QS_Tested = np.array([np.random.hypergeometric(Mistargeted_Quarantined[i], Population['Q_Recovered'][i,-1], I) if I >0 else 0 for i, I in enumerate(Anti_Cap)])
                Mistargeted_Quarantined = -QS_Tested + Mistargeted_Quarantined
                Anti_Cap -= QS_Tested
                QS_Exposed = np.random.binomial(QS_Tested, 1 - Anti_Test['Specificity'])

                QR_Tested = Anti_Cap
                QR_Exposed = np.random.binomial(QR_Tested, Anti_Test['Sensitivity'])
            else:
                QS_Exposed = 0
                QR_Exposed = 0
            I_Recovered = np.random.binomial(Population['Infectious'][:,-1] - I_Quarantined, gamma)
            QI_Recovered = np.random.binomial(Population['Q_Infectious'][:,-1], gamma)

            Population['Susceptible'] = np.hstack((Population['Susceptible'], np.array([Population['Susceptible'][:,-1] - S_Quarantined - S_Infected + QS_Exposed], dtype=int).T))
            Population['Infectious'] = np.hstack((Population['Infectious'], np.array([Population['Infectious'][:,-1] - I_Quarantined - I_Recovered + S_Infected], dtype=int).T))
            Population['Recovered'] = np.hstack((Population['Recovered'], np.array([Population['Recovered'][:,-1] + QR_Exposed + I_Recovered + QI_Recovered], dtype=int).T))
            Population['Q_Susceptible'] = np.hstack((Population['Q_Susceptible'], np.array([Population['Q_Susceptible'][:,-1] - QS_Exposed + S_Quarantined], dtype=int).T))
            Population['Q_Infectious'] = np.hstack((Population['Q_Infectious'], np.array([Population['Q_Infectious'][:,-1] - QI_Recovered + I_Quarantined], dtype=int).T))
            Population['Q_Recovered'] = np.hstack((Population['Q_Recovered'], np.array([Population['Q_Recovered'][:,-1] - QR_Exposed], dtype=int).T))

            if not Length is None:
                if t == Length:
                    break
            else:
                if all(np.sum([Population[K][:,-1] for K in ['Infectious','Q_Infectious']], 0) == 0):
                    break
            t += 1


        return Population