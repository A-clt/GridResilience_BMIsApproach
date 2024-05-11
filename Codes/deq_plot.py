import matplotlib.pyplot as plt
import pandas as pd
plt.rc('text',usetex=True)
import scienceplots
plt.style.use(['science','std-colors'])
plt.rcParams['figure.figsize'] = [0.99*7.01388889/2, 1.2]
plt.rcParams.update({'font.size': 6})

df1 = pd.read_csv('deq.csv')
mode1 = pd.read_csv('mode.csv')
df2 = pd.read_csv('deq2.csv')
mode2 = pd.read_csv('mode2.csv')
fig,ax = plt.subplots(2,1,sharex=True)
ax[0].plot(df1['Time'],df1['DEQ'],label = '$\\tau_a = 4$')
ax[0].plot(df2['Time'],df2['DEQ'],label = '$\\tau_a = 3$')
ax[0].set_ylabel('$\\|x-x_{eq}\\|^2$')
#list(df1['Time'])[-1]+0.01
ax[0].set_xlim([0,list(df1['Time'])[-1]+0.01])
#ax[0].set_ylim([0,3])
ax[0].legend(handlelength=1, frameon=True, fancybox=False, framealpha=0.8)
ax[0].grid()
ax[1].plot(mode1['Time'],mode1['MODE'])
ax[1].plot(mode2['Time'],mode2['MODE'])
ax[1].set_ylabel('Mode')
#axins.legend(handlelength=1, frameon=True, fancybox=False, framealpha=0.8)
ax[1].set_yticks([1,2,3,4,5,6])
ax[1].grid()  
ax[1].set_xlim([0,list(df1['Time'])[-1]+0.01])
ax[1].set_xlabel('Time [s]')
fig.align_ylabels()

plt.savefig('dist_eq_comp.pdf')



fig,ax = plt.subplots(1,1)
ax.plot(mode1['Time'],mode1['MODE'],label = '$\\tau_a = 4$')
ax.plot(mode2['Time'],mode2['MODE'],label = '$\\tau_a = 3$')
ax.set_ylabel('$\\|x-x_{eq}\\|^2$')
#list(df1['Time'])[-1]+0.01
ax.set_xlim([0,list(df1['Time'])[-1]+0.01])
ax.set_yticks([1,2,3,4,5,6])
ax.set_xlabel('Time [s]')
ax.legend(handlelength=1, frameon=True, fancybox=False, framealpha=0.8)
ax.grid()
plt.savefig('mode_comp.pdf')