import numpy as np
np.random.seed(8346)
import os
from pyomo.environ import *
import cvxpy as cv
import pyomo.environ as pyo
from dataclasses import dataclass
import matplotlib.pyplot as plt
plt.rc('text',usetex=True)
import scienceplots
plt.style.use(['science','std-colors'])
plt.rcParams['figure.figsize'] = [0.95*7.01388889, 2.2]
plt.rcParams.update({'font.size': 6})