#import numpy as np
from lib import *
from AttacksModel import Attacks
from constant import t_max,dt, t_cond,u_cond, eps1, eps2
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rc('text',usetex=True)
import scienceplots
plt.style.use(['science','std-colors'])
plt.rcParams['figure.figsize'] = [0.95*7.01388889, 2.2]
plt.rcParams.update({'font.size': 6})

class Simulator():
      def __init__(self):
            self.Attacks = Attacks()
            self.LTI = self.Attacks.LTI_dic[1]
            print(self.LTI.Hurwitz)
            print(np.linalg.eigvals(self.LTI.A))
            self.pstar = np.ones((self.LTI.ng,1))
            self.load = sum(self.pstar)
            self.build_u()
            self.xeq = np.linalg.inv(self.LTI.A)@(-self.LTI.B@self.u)
            self.x = self.xeq*(1+eps1)+eps2*np.ones(self.xeq.shape)
            self.x_history = list()
            self.xeq_history = list()
            self.u_history = list()
            self.mode_history = list()
            self.plot_mosaic = True

      def build_u(self):
            self.u = np.ones((self.LTI.ng+1,1))
            self.u[0] = self.load
            self.u[1:] = self.pstar

      def run(self):
            max_iter = int(t_max/dt)
            iter = 0
            cond_counter = 0
            init_load = self.load
            while iter<=max_iter:
                  self.Attacks.run()
                  #store results
                  self.x_history.append(self.x)
                  self.xeq_history.append(self.xeq)
                  self.u_history.append(self.u)
                  self.mode_history.append(self.Attacks.mode)
                  #update u
                  if abs(iter*dt-t_cond[cond_counter])<0.5*dt:
                        self.load = init_load + u_cond[cond_counter]
                        self.build_u()
                        cond_counter+= 1 if cond_counter<(len(t_cond)-1) else 0
                  self.x = self.Attacks.Ad@self.x + self.Attacks.Bd@self.u
                  self.xeq = np.linalg.inv(self.Attacks.A)@(-self.Attacks.B@self.u)                      
                  iter+=1      
      
      def plot(self):
            t = np.arange(0,t_max,dt)
            df = [self.x_history[i][0]*60 for i in range(int(t_max/dt))]
            load = [self.u_history[i][0] for i in range(int(t_max/dt))]
            pm = [[self.x_history[i][n] for i in range(int(t_max/dt))] for n in range(1,self.LTI.ng+1)]
            z = [self.x_history[i][-1] for i in range(int(t_max/dt))]
            deq = [np.linalg.norm(self.x_history[i]-self.xeq_history[i])**2 for i in range(int(t_max/dt))]
            mode = [self.mode_history[i] for i in range(int(t_max/dt))]
            print('Pc time in stable mode: %.3f'%(sum([1 if mode[i]==1 else 0 for i in range(len(mode))])/len(mode)*100))
            if not self.plot_mosaic:
                  _,ax = plt.subplots(1,1)
                  ax.plot(t,df, color='black')
                  ax.set_xlabel('Time [s]')
                  ax.set_ylabel('$\\Delta f$ [Hz]')
                  ax.grid()
                  ax.set_ylim([0.1*int(10*(min(df)-0.1)),0.2])
                  ax.set_xlim([0,t_max])
                  plt.savefig('Frequency.pdf')
                  
                  _,ax = plt.subplots(1,1)
                  for n in range(0,self.LTI.ng):
                        ax.plot(t,pm[n],label='$P^{m}_{%d}$'%(n+1))
                  ax.set_xlabel('Time [s]')
                  ax.set_ylabel('$P_m$ [pu]')
                  ax.grid()
                  ax.set_xlim([0,t_max])
                  plt.legend()      
                  plt.savefig('Power.pdf') 

                  _,ax = plt.subplots(1,1)
                  ax.plot(t,z, color='black')
                  ax.set_xlabel('Time [s]')
                  ax.set_ylabel('$z$ [pu]')
                  ax.grid()  
                  ax.set_xlim([0,t_max])  
                  plt.savefig('z.pdf') 

                  _,ax = plt.subplots(1,1)
                  ax.plot(t,deq, color='black')
                  ax.set_xlabel('Time [s]')
                  ax.set_ylabel('$\\|x-x_{eq}\\|^2$')
                  ax.grid()  
                  ax.set_xlim([0,t_max])  
                  plt.savefig('dist_eq.pdf')

                  _,ax = plt.subplots(1,1)
                  ax.plot(t,mode, color='black')
                  ax.set_xlabel('Time [s]')
                  ax.set_ylabel('Mode')
                  ax.set_yticks([1,2,3,4,5,6])
                  ax.grid()  
                  ax.set_xlim([0,t_max])  
                  plt.savefig('Mode.pdf')

            else:
                  
                  fig,ax_dict = plt.subplot_mosaic(
                  [
                        ["df", "pm"],
                        ["z", "deq"],
                        ["mode","mode"]
                  ],layout='constrained'
                  )   
                  ax_dict["df"].plot(t,df, color='black')
                  #ax_dict["df"].set_xlabel('Time [s]')
                  ax_dict["df"].set_ylabel('$\\Delta f$ [Hz]')
                  ax_dict["df"].grid()
                  ax_dict["df"].set_ylim([0.1*int(10*(min(df)-0.1)),0.1*int(10*(max(df)+0.1))])
                  ax_dict["df"].set_xlim([0,t_max])
                  
                  for n in range(0,self.LTI.ng):
                        ax_dict["pm"].plot(t,pm[n],label='$P^{m}_{%d}$'%(n+1), linestyle='solid' if n<7 else 'dotted')
                  #ax_dict["pm"].set_xlabel('Time [s]')
                  ax_dict["pm"].set_ylabel('$P_m$ [pu]')
                  ax_dict["pm"].grid()
                  ax_dict["pm"].set_xlim([0,t_max])
                  fig.legend(ncol=int(self.LTI.ng),loc='outside upper center',handlelength=1, frameon=True, fancybox=False, framealpha=0.8, fontsize=5)      

                  ax_dict["z"].plot(t,z, label='$z$', color='black')
                  ax_dict["z"].plot(t,load, label='$P_{\\mathrm{load}}$', color='red')
                  ax_dict["z"].set_xlabel('Time [s]')
                  ax_dict["z"].set_ylabel('$z/P_{\\mathrm{load}}$ [pu]')
                  ax_dict["z"].legend(ncol=2,handlelength=1, frameon=True, fancybox=False, framealpha=0.8)
                  ax_dict["z"].grid()  
                  ax_dict["z"].set_xlim([0,t_max])  
                  ax_dict["z"].sharex(ax_dict["df"])

                  ax_dict["deq"].plot(t,deq, color='black')
                  ax_dict["deq"].set_xlabel('Time [s]')
                  ax_dict["deq"].set_ylabel('$\\|x-x_{eq}\\|^2$')
                  ax_dict["deq"].grid()  
                  ax_dict["deq"].set_xlim([0,t_max]) 
                  ax_dict["deq"].sharex(ax_dict["pm"]) 
                  #df = pd.DataFrame({'Time':t,'DEQ':deq})
                  #df.to_csv('deq.csv', index=False)

                  """
                  axins = ax_dict["df"].inset_axes([0.11, 0.13, 0.85, 0.55])
                  axins.plot(t,mode, label='$\\tau_{a}=%.1f$'%(self.Attacks.tau_a), color='black')
                  #axins.set_xlabel('Time [s]')
                  axins.set_ylabel('Mode')
                  #axins.legend(handlelength=1, frameon=True, fancybox=False, framealpha=0.8)
                  axins.set_yticks([1,2,3,4,5,6])
                  axins.grid()  
                  axins.set_xlim([0,t_max])

                  """
                  ax_dict["mode"].plot(t,mode, label='$\\tau_{a}=%.1f$'%(self.Attacks.tau_a), color='black')
                  ax_dict["mode"].set_xlabel('Time [s]')
                  ax_dict["mode"].set_ylabel('Mode')
                  ax_dict["mode"].set_yticks([1,2,3,4,5,6])
                  ax_dict["mode"].grid()  
                  ax_dict["mode"].set_xlim([0,t_max])  
                  df = pd.DataFrame({'Time':t,'MODE':mode})
                  df.to_csv('mode.csv', index=False)
                  # cons_mode = list()
                  # cons_mode.append(mode[0])
                  # width = list()
                  # start_t = list()
                  # start_t.append(0)
                  # m_previous = mode[0]
                  # width_counter = 0
                  # for m in mode[1:]:
                  #       if m!=m_previous:
                  #             m_previous=m
                  #             cons_mode.append(m)
                  #             width.append(width_counter*dt)
                  #             start_t.append(start_t[-1]+width[-1])
                  #             width_counter=0
                  #       else:
                  #             width_counter+=1      
                  # width.append(t[-1]-sum(width))
                  # print(cons_mode)
                  # print(width)
                  # for i,m in enumerate(cons_mode):
                  #       # Make the rectangle
                  #       color_rect = 'green' if m==1 else 'red'
                  #       rect = Rectangle(xy=(start_t[i], 0),
                  #                               width=width[i],
                  #                               height=1,
                  #                               color=color_rect)

                  #       # Add Rectangle to the figure
                  #       ax_dict["mode"].add_patch(rect)
                  #       #ax_dict["mode"].plot(t,mode, label='$\\tau_{a}=%.1f$'%(self.Attacks.tau_a), color='black')
                  # ax_dict["mode"].set_xlabel('Time [s]')
                  # ax_dict["mode"].set_ylabel('Mode')
                  # #ax_dict["mode"].set_yticks([1,2,3,4,5,6])
                  # ax_dict["mode"].tick_params(left=False,right=False,labelleft=False,which='both')
                  # #ax_dict["mode"].grid()  
                  # ax_dict["mode"].set_xlim([0,t_max])
                  #fig.legend(loc='outside upper center',handlelength=1, frameon=True, fancybox=False, framealpha=0.8, fontsize=5)      


                  fig.align_ylabels()
                  plt.savefig('Mosaic_Stable_VaryingLoad_v2.pdf')                  


if __name__=='__main__':
      sim = Simulator()
      sim.run()
      sim.plot()