import numpy as np
dt = 0.001
t_max = 500
#Time varying Load
t_cond = np.arange(0,t_max,0.01)
u_cond = 5*np.sin(0.01*t_cond)
eps1 = 0
eps2 = 0
##Perturb IC
#t_cond = [0]
#u_cond = [0]
#eps1 = 0.05
#eps2 = 0.0005
##Transient
#t_cond = [10]
#u_cond = [5]
#eps1 = 0
#eps2 = 0