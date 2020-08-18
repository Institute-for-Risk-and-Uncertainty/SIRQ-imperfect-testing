# SIRQ-Model
Binomial model for SIR model with quarantine

# Recreating the figures from the paper
The scratch script can be run to recreate figures 5-8 from the paper, saving each as an eps file. Creation of each can be disabled by setting the flags to False.

# Running a new scenario
Model parameters must be initialised for both the antigen and antibody tests. The model must first be created as an object using the command:
``` python
QTest = SIRQModel.SIRQ()
```

Antibody and antigen tests can then be defined using the following commands:

``` python
QTest.Def_Anti_Test(0.9, 0.9, 1e5, 0.8, 1)
QTest.Def_Inf_Test(0.9, 0.9, 1e5, 0.8, 1)
```

This initialises both tests with (1) sensitivity = 0.9, (2) specificity = 0.9, (3) capacity = 1e5, (4) targeting = 0.8, and (5) test interval = 1.

The target population must then be initialised:

``` python
PopSize = 6.7e7
QTest.Def_Population(0.984*PopSize, 0.01*PopSize, 0.001*PopSize, 0*PopSize, 0.004*PopSize, 0.001*PopSize)
```

This initialises the population with a split of: S = 98.4%, I = 1%, R = 0.1%, Q_S = 0%, Q_I = 0.4%, and Q_R = 0.1% of the population.

A model run can then be performed with SIR parameters beta and gamma:

``` python
Pop = QTest.Run(0.32, 0.1, Length = 100, Runs = 1)
```

This performs a single model evaluation with beta=0.32 and gamma=0.1 for 100 days. This returns a dictionary, 'Pop', which contains a list of integers describing the population within each state for each time step. Each state list can be accessed using the key for that state: 'Susceptible', 'Infectious', 'Recovered', 'Q_Susceptible', 'Q_Infectious', and 'Q_Recovered'.
