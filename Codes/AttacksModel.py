#import numpy as np
from lib import *
from BuildSystems import LTI_system
from constant import dt
tau_z = 10
class Attacks():
      def __init__(self):
            self.tau_d = 2
            self.N0 = 5
            self.tau_a = 4
            self.T0 = 2
            self.LTI_dic = {
                  1: LTI_system(40,1.5,tau_z,alpha='EquivalentSharing'),
                  2: LTI_system(40,-100,tau_z,alpha='EquivalentSharing'),
                  3: LTI_system(40,-200,tau_z,alpha='EquivalentSharing'),
                  4: LTI_system(40,-150,tau_z,alpha='EquivalentSharing'),
                  5: LTI_system(40,-120,tau_z,alpha='EquivalentSharing'),                  
                  6: LTI_system(40,-170,tau_z,alpha='EquivalentSharing'),
            }
            self.Q = [1,2,3,4,5,6]
            self.Qu = [2,3,4,5,6]
            self.mode = 1
            self.adttimer = 0
            self.aattimer = 0
            self.flow = 0
            self.jump = 0
            self.Ad = self.LTI_dic[1].Ad
            self.Bd = self.LTI_dic[1].Bd  
            self.A = self.LTI_dic[1].A
            self.B = self.LTI_dic[1].B  
            for i in range(1,7):
                  np.savetxt('A%d.txt'%i,self.LTI_dic[i].A)             
            np.savetxt("B.txt",self.B)
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
