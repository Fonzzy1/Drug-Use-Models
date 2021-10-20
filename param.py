import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def SUAR_norm(p1,p2,p3, p4, p5,d,S0,U0,A0,R0,i,h):

    vect_initial = np.array([S0,U0,A0,R0]).T

    j = 0
    vect = vect_initial.copy()

    ds = d[0]
    du = d[1]
    da = d[2]
    dr = d[3]

    while j < i:
        vectplus1 = np.zeros([4])
        s = vect[0]
        u = vect[1]
        a = vect[2]
        r = vect[3]
        vectplus1[0] = s +h*((-s*p1*(u+a)+u*p2*(s+r)  -ds*s + 1)-s *(1 - ds*s - du*u - da*a - dr*r))
        vectplus1[1] = u +h*((s*p1*(u+a)-u*p2*(s+r) - u*p3 -du*u )-u*(1 - ds*s - du*u - da*a - dr*r))
        vectplus1[2] = a +h*((u*p3 -a*p4*(s+r) + r*p5*(u+a) -da*a )-a*(1 - ds*s - du*u - da*a - dr*r))
        vectplus1[3] = r +h*((a*p4*(s+r) - r*p5*(u+a) -dr*r )-r*(1 - ds*s -du*u - da*a - dr*r))
        vect = vectplus1.copy()
        j+=1

    if s + r > 0.98:
        return True
    else:
        return False

res = pd.DataFrame(columns=['p1','p2','p3','p4','p5', 'extinct'])
for p1 in range(-3,4):
    for p2  in range(-3,4):
        for p3 in range(-3,4):
            for p4 in range(-3,4):
                for p5 in range(-3,4):

                    bool = SUAR_norm(2**p1, 2**p2, 2**p3, 2**p4, 2**p5, [0.1, 0.5, 1, 0.1], 0.6, 0.2, 0.1, 0.1, 50, 0.1)
                    res.loc[len(res)] = [float(p1),float(p2),float(p3),float(p4), float(p5), bool]

p = ['p1','p2','p3','p4','p5']

xlabels = ['-3','-2','-1','0','1','2','3']


fig, axs = plt.subplots(5, 5)
for i , a in enumerate( p):
    for j , b in enumerate(p):
        if i !=j:
            w = res[[a,b,'extinct']].groupby([a,b]).sum().unstack()
            t = axs[i,j].pcolormesh(w, vmin = 0, vmax = 216)
        if i == 4:
            axs[i,j].set_xticks(np.arange(7) + 0.5, minor=False)
            axs[i,j].set_xticklabels(xlabels, fontdict=None, minor=False)
            axs[i,j].set_xlabel('Log2(' + b +')')
        else:
            axs[i,j].get_xaxis().set_visible(False)
        if j == 0 :
            axs[i,j].set_yticks(np.arange(7) + 0.5, minor=False)
            axs[i,j].set_yticklabels(xlabels, fontdict=None, minor=False)
            axs[i,j].set_ylabel('Log2(' + a +')')
        else:
            axs[i,j].get_yaxis().set_visible(False)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(t, cax=cbar_ax)
fig.suptitle("Number of Simulations that reach extinction")
plt.show()
