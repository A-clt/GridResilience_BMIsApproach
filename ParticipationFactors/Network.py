#import numpy as np
from lib import *
from constant import Kg, Tsm

Sbase = 100
Vbase = 230
omega_syn = 2*np.pi*60
nb = 6
ng = 4

@dataclass
class Generator:
    bus: int
    Sgen: float
    H: float
    D: float
    Tsm: float
    Kg: float

@dataclass
class Line:
    bus_from: int
    bus_to: int
    R: float
    X: float 
    B: float 
    gen_from: int
    gen_to: int  

@dataclass
class Transformer:
    bus_from: int
    bus_to: int
    R: float
    X: float 
    gen_from: int
    gen_to: int

def create_dic_generator():
      dic_gen = {'Gen1':Generator(bus=1,Sgen=1000,H=5,D=1.5, Tsm=Tsm, Kg=Kg),
                 'Gen2':Generator(bus=2,Sgen=1000,H=5,D=1.5, Tsm=Tsm, Kg=Kg),
                 'Gen3':Generator(bus=3,Sgen=1000,H=5,D=1.5, Tsm=Tsm, Kg=Kg),
                 'Gen4':Generator(bus=6,Sgen=1000,H=5,D=1.5, Tsm=Tsm, Kg=Kg)
                 }
      return dic_gen

def create_dic_line():
      dic_line = {'Line1':Line(bus_from=1,bus_to=2,R=0.0035,X=0.0411,B=0.6987, gen_from=1, gen_to=2),
                  'Line2':Line(bus_from=2,bus_to=3,R=0.001,X=0.025,B=0.75, gen_from=2, gen_to=3),
                  'Line3':Line(bus_from=2,bus_to=4,R=0.0013,X=0.0151,B=0.2572, gen_from=2, gen_to=4),
                  'Line4':Line(bus_from=4,bus_to=5,R=0.007,X=0.0086,B=0.146, gen_from=4, gen_to=5),
                  'Line5':Line(bus_from=4,bus_to=6,R=0.0013,X=0.0213,B=0.2214, gen_from=4, gen_to=6)
                 }   
      return dic_line  

def create_dic_transformer():
      dic_transformer = {}
      return dic_transformer

def load6bus():
     gen = create_dic_generator()
     line = create_dic_line()
     trafos = create_dic_transformer()
     return gen,line,trafos,Sbase,ng,nb
