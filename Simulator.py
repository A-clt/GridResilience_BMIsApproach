from lib import *
from AttacksModel import Attacks
from constant import t_max,dt, t_cond,u_cond, eps1, eps2, filename

class Simulator():
      def __init__(self):
            self.Attacks = Attacks()
            self.LTI = self.Attacks.LTI_dic[1]
            self.pstar = np.ones((self.LTI.ng,1))
            self.load = sum(self.pstar)
            self.build_u()
            self.xeq = np.linalg.inv(self.LTI.A)@(-self.LTI.B@self.u)
            self.x = self.xeq*(1+eps1)+eps2*np.ones(self.xeq.shape)
            self.x_history = list()
            self.xeq_history = list()
            self.u_history = list()
            self.mode_history = list()

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
          
            fig,ax_dict = plt.subplot_mosaic(
            [
                  ["df", "pm"],
                  ["z", "deq"],
                  ["mode","mode"]
            ],layout='constrained'
            )   
            ax_dict["df"].plot(t,df, color='black')
            ax_dict["df"].set_ylabel('$\\Delta f$ [Hz]')
            ax_dict["df"].grid()
            ax_dict["df"].set_ylim([0.1*int(10*(min(df)-0.1)),0.1*int(10*(max(df)+0.1))])
            ax_dict["df"].set_xlim([0,t_max])
            
            for n in range(0,self.LTI.ng):
                  ax_dict["pm"].plot(t,pm[n],label='$P^{m}_{%d}$'%(n+1), linestyle='solid' if n<7 else 'dotted')
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

            ax_dict["mode"].plot(t,mode, label='$\\tau_{a}=%.1f$'%(self.Attacks.tau_a), color='black')
            ax_dict["mode"].set_xlabel('Time [s]')
            ax_dict["mode"].set_ylabel('Mode')
            ax_dict["mode"].set_yticks([1,2,3,4,5,6])
            ax_dict["mode"].grid()  
            ax_dict["mode"].set_xlim([0,t_max])  

            fig.align_ylabels()
            if not os.path.exists('Figures'):
                  os.makedirs('Figures')            
            plt.savefig('Figures//%s.pdf'%filename)                  

if __name__=='__main__':
      sim = Simulator()
      sim.run()
      sim.plot()