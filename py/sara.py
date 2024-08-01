"""
Usage:
	pretty much have to hard-code in K,s, and the order of the state variables
	Call function `sara` at each time step, passing in the current sim time and the list of new state space variables it needs
	That in turn calls the specific solver method for each equation
"""
from math import exp
class Solver:
	prev_times:list[float]
class L:
	K:list
	s:list
	offset:list
	def q1(self,i:int):
		q0=(delta_n/zeta(i))*(1-Phi(i))
		return q0,
	def Phi(self,i,delta):
		return exp(self.s[i]*delta)
	def zeta(self,i,delta):
		return -(self.s[i]*delta)
	def __init__(self,K,s,offset):
		self.K=K
		self.s=s
		self.offset=offset
	def __call__(self,solv:Solver):

def calc_L(solv:Solver):
	K=[1.8e-9,complex(2.703e-16,7.316e-16), complex(2.703e-16,7.316e-16)]
	s=[complex(7.108e5,-9.464e5), complex(7.108e5,9.464e5)]
	scale=1.8e-9
	scale_offset=0
	temp:float = 0
	for i in range(len(K)):
		if K[i] == K[i].conjugate():



T=40
delta_n=1
y=[]
x=[]
R2=2

def q2(i:int):
	zi=zeta(i)
	q0=(delta_n/zi**2)*(-1+zi+Phi(i))
	q1=(delta_n/zi**2)*(1-(1+zi)*Phi(i))
	return q0,q1
	prev_state_vars:list[list[float]]
def sara(t_new:float,x:list[float]):
