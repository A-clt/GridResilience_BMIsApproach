#import numpy as np
from lib import *
from Network import loadIEEE39bus
from constant import dt, beta

class LTI_system():
      def __init__(self,M_DERS,D_DERS,Tz,alpha='ProportionalSharing'):
            self.gen,self.line,self.trafos,self.Sbase,self.ng,self.nb = loadIEEE39bus()
            self.M_DERS = M_DERS
            self.D_DERS = D_DERS
            self.Tz = Tz
            self.set_alpha(alpha)
            self.build_params_mat()
            self.beta = beta
            self.build_A_mat()
            self.build_B_mat()
            self.c2d()

      def set_alpha(self,alpha):
            if alpha=='EquivalentSharing':
                  self.alpha = 1/self.ng*np.ones((self.ng,1))
            elif alpha=='ProportionalSharing':
                  self.alpha = np.ones((self.ng,1))
                  i = 0
                  for _,value in self.gen.items():
                        self.alpha[i] = value.Sgen
                  self.alpha = self.alpha/sum(self.alpha)
            elif alpha=='RandomDispatch':
                  self.alpha = np.array([[0.05],[0.1],[0.15],[0.08],[0.12],[0.3],[0.05],[0.025],[0.075],[0.05]]) 
                  self.alpha = self.alpha/sum(self.alpha)
                  print(self.alpha)
                  input()     
            else:
                  self.alpha=alpha      

      def build_params_mat(self):
            M = self.M_DERS*np.ones((self.nb,1))
            D = self.D_DERS*np.ones((self.nb,1))
            self.Tg = np.zeros((self.ng,self.ng))
            self.Kp = np.zeros((self.ng,1))
            idx=0
            for _,value in self.gen.items():
                  M[value.bus-1] = 2*value.H*value.Sgen/self.Sbase
                  D[value.bus-1] = value.D
                  self.Tg[idx][idx] = value.Tsm
                  self.Kp[idx] = value.Kg*value.Sgen/self.Sbase 
                  idx+=1
            self.M_eff = sum(M)
            self.D_net = sum(D)   

      def build_A_mat(self):
            self.A = np.block([
                  [-self.D_net/self.M_eff,1/self.M_eff*np.ones((1,self.ng)),0],
                  [-np.linalg.inv(self.Tg)@self.Kp,-np.linalg.inv(self.Tg),np.linalg.inv(self.Tg)@self.alpha],
                  [1/self.Tz*self.beta,1/self.Tz*np.ones((1,self.ng)),-1/self.Tz]
            ])

            self.eig,_ = np.linalg.eig(self.A)
            if np.array([e<0 for e in np.real(self.eig)]).all():
                  self.Hurwitz = True
            else:
                  self.Hurwitz = False                

      def build_B_mat(self): 
            self.B = np.block([[np.block([-1/self.M_eff,np.zeros((1,self.ng))])],
                        [np.block([np.zeros((self.ng,1)),np.linalg.inv(self.Tg)@(np.eye(self.ng)-self.alpha*np.ones((1,self.ng)))])],
                        [np.block([np.zeros((1,1+self.ng))])]])     
         
      def c2d(self):
            n = self.A.shape[0]  # Assuming square matrices

            M1 = np.eye(n) + 0.5*dt * self.A
            M2 = np.eye(n) - 0.5*dt * self.A

            self.Ad = M1@np.linalg.inv(M2)
            self.Bd = 0.5*dt*(np.eye(n)+self.Ad)@self.B