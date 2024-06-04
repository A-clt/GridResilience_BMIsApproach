from lib import *
from BuildSystems import LTI_system
from constant import dt, tau_z, Part_Factors, tau_d, N0, tau_a, T0, M_DERS, D_DERS

class Attacks():
      def __init__(self):
            self.tau_d = tau_d
            self.N0 = N0
            self.tau_a = tau_a
            self.T0 = T0
            self.LTI_dic = {
                  1: LTI_system(M_DERS[0],D_DERS[0],tau_z,alpha=np.array([[0.2],[0.4],[0.2],[0.2]])),
                  2: LTI_system(M_DERS[1],D_DERS[1],tau_z,alpha=np.array([[0.5], [1e-10], [0.5], [1e-10]])),
                  3: LTI_system(M_DERS[2],D_DERS[2],tau_z,alpha=np.array([[0.33], [0.33], [1e-10], [0.33]]))
            }
            self.Q = [1,2,3]
            self.Qu = [2,3]
            self.mode = 1
            self.adttimer = 0
            self.aattimer = 0
            self.flow = 0
            self.jump = 0
            self.Ad = self.LTI_dic[1].Ad
            self.Bd = self.LTI_dic[1].Bd  
            self.A = self.LTI_dic[1].A
            self.B = self.LTI_dic[1].B  
            if not os.path.exists('LTIs'):
                  os.makedirs('LTIs')
            for i in range(1,4):
                  np.savetxt('LTIs//A%d.txt'%i,self.LTI_dic[i].A)             
            np.savetxt("LTIs//B.txt",self.B)
            
      def C(self):
            if (0<=self.adttimer and self.adttimer<=self.N0) and (0<=self.aattimer and self.aattimer<=self.T0):
                  self.flow=1
            else:
                  self.flow=0

      def D(self):
            if (1<=self.adttimer and self.adttimer<=self.N0) and (0<=self.aattimer and self.aattimer<=self.T0):
                  self.jump=1
            else:
                  self.jump=0

      def G(self):
            if self.mode==1:
                  self.mode=np.random.choice(np.array(self.Qu)) 
            else:
                  self.mode=1      
            self.adttimer = max(self.adttimer-1,0)

      def run(self):
            self.C()
            self.D()
            if self.flow:
                  self.Ad = self.LTI_dic[self.mode].Ad
                  self.Bd = self.LTI_dic[self.mode].Bd
                  self.A = self.LTI_dic[self.mode].A
                  self.B = self.LTI_dic[self.mode].B                  
                  self.ind = 1 if self.mode in self.Qu else 0
                  temp_adttimer = self.adttimer 
                  temp_aattimer = self.aattimer 
                  self.adttimer+= dt*1/self.tau_d
                  self.aattimer+= dt*(1/self.tau_a-self.ind)
                  self.C()
                  if not self.flow:
                        self.jump = 1
                        self.adttimer = temp_adttimer
                        self.aattimer = temp_aattimer     
            if self.jump:
                  self.G()    
