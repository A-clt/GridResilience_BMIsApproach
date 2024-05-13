from AttacksModel import Attacks
from lib import *

# for pyomo
def array_to_dict(array):
    rows, cols = array.shape
    result_dict = {}
    for i in range(rows):
        for j in range(cols):
            result_dict[(i, j)] = array[i, j]
    return result_dict

class LMI_solver():
      def __init__(self):
            self.At = Attacks()
            self.Au1 = self.At.LTI_dic[2].A
            self.Au2 = self.At.LTI_dic[3].A
            self.Au3 = self.At.LTI_dic[4].A
            self.Au4 = self.At.LTI_dic[5].A
            self.Au5 = self.At.LTI_dic[6].A
            self.A = self.At.LTI_dic[1].A
            self.m = pyo.ConcreteModel()
            self.tau_a = list()
            self.c_prime_list = list()
            self.c_list = list()
            self.solve_model_CVXPY()

      def model_CVXPY(self):
            # create optimization variable
            self.P = cv.Variable((self.A.shape[0], self.A.shape[1]))
            self.c = cv.Parameter(nonneg=True)
            self.c_prime = cv.Parameter(nonneg=True)

            # create constraints
            PD_cstr = self.P-1e-6*np.eye(self.A.shape[0])>>0
            P_sym_cstr = self.P==self.P.T
            P_stab_cstr = self.A.T@self.P+self.P@self.A+self.c*self.P<<0
            P_ustab1_cstr = self.Au1.T@self.P+self.P@self.Au1-self.c_prime*self.P<<0
            P_ustab2_cstr = self.Au2.T@self.P+self.P@self.Au2-self.c_prime*self.P<<0
            P_ustab3_cstr = self.Au3.T@self.P+self.P@self.Au3-self.c_prime*self.P<<0
            P_ustab4_cstr = self.Au4.T@self.P+self.P@self.Au4-self.c_prime*self.P<<0
            P_ustab5_cstr = self.Au5.T@self.P+self.P@self.Au5-self.c_prime*self.P<<0

            # create optimization problem
            self.optprob = cv.Problem(cv.Minimize(0),constraints=[PD_cstr,P_sym_cstr,P_stab_cstr,P_ustab1_cstr,P_ustab2_cstr,P_ustab3_cstr,P_ustab4_cstr,P_ustab5_cstr])
      
      def solve(self):
            try:
                  self.optprob.solve('CLARABEL')
                  if self.optprob.status=='optimal':
                        if self.check_inequalities(self.P.value,self.c.value,self.c_prime.value):
                              self.c_prime_list.append(self.c_prime.value)
                              self.c_list.append(self.c.value)
                              tau_a = (self.c_prime.value+self.c.value)/(self.c.value)
                              self.tau_a.append(tau_a) 
            except cv.error.SolverError:
                  pass                     

      def solve_model_CVXPY(self):
            c = np.arange(0.1,1,0.1)
            c_prime = np.arange(1,10,0.25)
            self.model_CVXPY()
            for j in range(len(c)):
                  for i in range(len(c_prime)):
                        if (1+c_prime[i]/c[j])>100:
                              pass
                        else:
                              self.c.value = c[j]
                              self.c_prime.value = c_prime[i]
                              self.solve()                
            minpos = self.tau_a.index(min(self.tau_a))      
            print('Bound on tau_a:%.3f, value of lambda_s:%.3f, value of lambda_u:%.3f'%(min(self.tau_a),self.c_list[minpos],self.c_prime_list[minpos]))

      def check_inequalities(self,P,c,c_prime):
            cd1 = np.all(np.linalg.eigvals(P) > 0)  
            cd2 = np.all(np.linalg.eigvals(self.A.T@P + P@self.A + c*P)<=0)
            cd3 = np.all(np.linalg.eigvals(self.Au1.T@P + P@self.Au1 - c_prime*P)<=0)
            cd4 = np.all(np.linalg.eigvals(self.Au2.T@P + P@self.Au2 - c_prime*P)<=0)
            cd5 = np.all(np.linalg.eigvals(self.Au3.T@P + P@self.Au3 - c_prime*P)<=0)
            cd6 = np.all(np.linalg.eigvals(self.Au4.T@P + P@self.Au4 - c_prime*P)<=0)
            cd7 = np.all(np.linalg.eigvals(self.Au5.T@P + P@self.Au5 - c_prime*P)<=0)
            return cd1 and cd2 and cd3 and cd4 and cd5 and cd6 and cd7

if __name__=='__main__':
     LMI = LMI_solver()

                         