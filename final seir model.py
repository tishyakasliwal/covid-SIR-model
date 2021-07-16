



from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import turicreate



def deriv(y, t, N, beta, gamma, delta, alpha_opt, rho):
    S, E, I, R, D = y
    
    def alpha(t):
       
        return s * I/N + alpha_opt
    
    dSdt = -beta(t) * S * I / N
    dEdt = beta(t) * S * I / N - delta * E
    dIdt = delta * E - (1 - alpha(t)) * gamma * I - alpha(t) * rho * I
    
    dRdt = (1 - alpha(t)) * gamma * I
    dDdt = alpha(t) * rho * I
    return dSdt, dEdt, dIdt, dRdt, dDdt


L=55
Second=75
Third= 95
Fourth = 109
Fifth = 122
Sixth = 133
Seventh = 165
N = 1_387_297_452 
D = 14 # infections lasts four days
gamma = 1/D
delta = 1.0 / 5  # incubation period of five days
s = 0.02
alpha = 0.025
def R_0(t):
    if t<L:
        return 3.5
    elif t < Second and t>L:
        return 1.8# 1.71
    elif t<Third and t>Second:
        return 1.5 #1.46
    elif t<Fourth and t>Third:
        return 1.3 #1.27
    elif t<Fifth and t>Fourth:
        return 1.25#1.23
    elif t<Sixth and t>Fifth:
        return 1.22
    elif t> Sixth and t<Seventh:
        return 1.22
    elif t> Seventh:
        return 3
def beta(t):
    return R_0(t) * gamma

alpha_opt = 0.035  # 5% death rate
rho = 1/18  # 9 days from infection until death
S0, E0, I0, R0, D0 = N-1, 1, 0, 0, 0  # initial conditions: one exposed


# In[22]:


t = np.linspace(0, 499, 500) # Grid of time points (in days)
y0 = S0, E0, I0, R0, D0 # Initial conditions vector

# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta, alpha_opt, rho))
S, E, I, R, D = ret.T


# In[23]:




sf = turicreate.SFrame("covid india final.csv")




fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(sf['Date'], sf['Total cases'], 'b', alpha=0.5, lw=2, label='Infected actual')
ax.plot(sf['Date'], sf['Deaths'], 'g', alpha=0.5, lw=2, label='Deaths actual')
ax.plot(t, R, 'g', alpha=0.5, lw=2, label='Recovered')
ax.plot(t, E, 'y', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, D, 'k', alpha=0.5, lw=2, label='Dead')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number(thousands)')

ax.set_ylim(0,1000000)

ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()


import plotly.graph_objects as go
import numpy as np

t = np.linspace(0, 499, 500)
# Create traces
fig = go.Figure()
fig.add_trace(go.Scatter(x=t, y=S,
                    mode='lines',
                    name='susceptible'))
fig.add_trace(go.Scatter(x=t, y=I, mode='lines+markers', name='infected'))
fig.add_trace(go.Scatter(x=t, y=R, mode='markers', name='recovered'))
fig.add_trace(go.Scatter(x=t, y=E, mode='lines', name='exposed'))
fig.add_trace(go.Scatter(x=t, y=D,
                    mode='lines',
                    name='dead'))

fig.show()





import plotly.express as px




import plotly.offline
import plotly.graph_objs as go
plotly.offline.plot(fig, filename='seird.html')





fig.write_html("seird.html")

