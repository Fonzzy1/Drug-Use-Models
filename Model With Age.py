import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

C = pd.DataFrame(columns=['s', 't', 'a', 'r'])
T = pd.DataFrame(columns=['s', 't', 'a', 'r'])
Y = pd.DataFrame(columns=['s', 't', 'a', 'r'])
A = pd.DataFrame(columns=['s', 't', 'a', 'r'])
S = pd.DataFrame(columns=['s', 't', 'a', 'r'])
"""
Society Factors

Aging rate is the portion of people who advance in age each step
Birthrate is determined by sum number of adults of birthing age
deathrate the proportion of that group that dies each iteration
"""
aging_rate = 0.1
birth_rate = 0.09
death_rate = 0.1*np.array([[0.2, 0.2, 0.2, 0.2],
                       [0.05, 0.05, 0.05, 0.05],
                       [0.1, 0.1, 0.1, 0.1],
                       [0.2, 0.2, 0.2, 0.2],
                       [0.5, 0.5, 0.5, 0.5]])

"""
X0 is the matrix that contains the initial data:

[[CS, CT, CA, CR],
 [TS, TT, TA, TR],
 [YS, YT, YA, YR],
 [AS, AT, AA, AR],
 [SS, ST, SA, SR]]
"""

X0 = 100 * np.array([[5, 0, 0, 0],
                     [10, 10, 0, 0],
                     [10, 5, 5, 0],
                     [10, 0, 100, 0],
                     [10, 0, 0, 0]])

(n_age_groups, n_status) = np.shape(X0)

"""
Interaction Matrix - How much each age group interacts/influences with other age groups
Doesnt necessarily need to be symentrical eg, adults influence kids more than kids influence adults
I[i,j] = effect on age group i from group j


Param Matrix - This is a nice way to store the parameters for each age_group

This is where out key assumptions come in
Going across the row:

Movement from S to T, based on interactions with A+T and age influence
Movement from T to S, based on interaction with S+R and age influence
Movement from T to A
Movement from R to A, based on interactions with A+T and age influence
Movement from A to R, based on interaction with S+R and age influence
"""

I =  0.001*np.array([[0.1, 0.0, 0.2, 0.5, 0.2],
                       [0.0, 0.5, 0.2, 0.2, 0.1],
                       [0.0, 0.1, 0.7, 0.2, 0.0],
                       [0.0, 0.0, 0.1, 0.8, 0.1],
                       [0.0, 0.0, 0.0, 0.5, 0.5]])

P = 0.1*np.array([[0.0, 1.0, 0.9, 0.0, 1.0],
              [0.5, 0.5, 0.9, 0.8, 0.2],
              [0.5, 0.5, 0.9, 0.6, 0.2],
              [0.1, 0.1, 0.9, 0.6, 0.1],
              [0.0, 0.0, 0.9, 0.9, 0.1]])




"""
Now time to simulate the population
This done in two distinct parts, first implementing the change between groups sideways and then adding the aging effect
"""

itterations  = 0
while itterations < 100:
    C.loc[itterations] = X0[0]
    T.loc[itterations] = X0[1]
    Y.loc[itterations] = X0[2]
    A.loc[itterations] = X0[3]
    S.loc[itterations] = X0[4]

    i = 0

    IE = np.matmul(I, X0)

    Xn = np.zeros(np.shape(X0))

    while i < n_age_groups:
        group_vect = X0[i].copy()
        s = group_vect[0]
        t = group_vect[1]
        a = group_vect[2]
        r = group_vect[3]
        interaction_neg = np.sum(IE[i][[1, 2]])
        interaction_pos = np.sum(IE[i][[0, 3]])

        Xn[i][0] = s - (s * P[i][0] * interaction_neg) + t * P[i][1] * interaction_pos
        Xn[i][1] = t + (s * P[i][0] * interaction_neg) - t * P[i][1] * interaction_pos - t * P[i][2]
        Xn[i][2] = a + t * P[i][2] + r * P[i][3] * interaction_neg - a * P[i][4] * interaction_pos
        Xn[i][3] = r - r * P[i][3] * interaction_neg + a * P[i][4] * interaction_pos

        i += 1

    Xn_aged = np.zeros(np.shape(X0))

    i = 0
    j = 0

    while i < n_age_groups:
        j = 0
        while j < n_status:
            if i == 0:
                if j == 0:
                    Xn_aged[i][j] = Xn[i][j] * (1 - aging_rate - death_rate[i][j]) + birth_rate * sum(Xn[2] + Xn[3])
                else:
                    Xn_aged[i][j] = Xn[i][j] * (1 - aging_rate - death_rate[i][j])
            else:
                Xn_aged[i][j] = Xn[i][j] * (1 - aging_rate - death_rate[i][j]) + X0[i - 1][j] * aging_rate
            j += 1
        i += 1

    X0 = Xn_aged
    itterations += 1


plt.plot(A.index.values,C['s'])
plt.plot(A.index.values,C['t'])
plt.plot(A.index.values,C['a'])
plt.plot(A.index.values,C['r'])
plt.plot(A.index.values,T['s'])
plt.plot(A.index.values,T['t'])
plt.plot(A.index.values,T['a'])
plt.plot(A.index.values,T['r'])
plt.plot(A.index.values,Y['s'])
plt.plot(A.index.values,Y['t'])
plt.plot(A.index.values,Y['a'])
plt.plot(A.index.values,Y['r'])
plt.plot(A.index.values,A['s'])
plt.plot(A.index.values,A['t'])
plt.plot(A.index.values,A['a'])
plt.plot(A.index.values,S['r'])
plt.plot(A.index.values,S['s'])
plt.plot(A.index.values,S['t'])
plt.plot(A.index.values,S['a'])
plt.plot(A.index.values,S['r'])
plt.show()

