"""
Usage:
	pretty much have to hard-code in K,s, and the order of the state variables
	Call function `sara` at each time step, passing in the current sim time and the list of new state space variables it needs
	That in turn calls the specific solver method for each equation
"""
from math import exp, sin
from pade import Pade,separate
from poly import Poly
from matplotlib import pyplot as plt
from random import uniform
from typing import Callable
import numpy as np

class Solver:
	prev_time:float
	prev_delta_ts:list[float]
	prev_values:list[list[float]]
	num_iter:int
	K:list
	s:list
	idxs:list[str]
	q:Callable
	def __init__(self,K:list,s:list):
		self.prev_time=0
		self.num_iter=0
		self.K=K
		self.s=s
		self.prev_delta_ts=[0]*8
		self.prev_values=[[0]*5 for i in range(len(K)+2)]
		self.idxs=['y'+str(i) for i in range(len(K))] + ['y','x']
		if len(K)==1:
			self.q=q1
		elif len(K)==2:
			self.q=q2
		elif len(K)==3:
			self.q=q3
		else:
			self.q=q4
	def __repr__(self):
		s=f"""{self.num_iter} iterations, final state variables are:
			{self.prev_values[self.idxs.index('y')][-1]}
			{len(self.prev_values[0])} y0's
			{len(self.prev_values[1])} y1's
			"""
		return s
	def update(self,t,x):
		x_idx = self.idxs.index('x')
		prev_x = self.prev_values[x_idx][0]
		delta_t = t-self.prev_time
		new_x = prev_x + (delta_t*x)
		self.prev_values[x_idx] = [new_x] + self.prev_values[x_idx]
		self.prev_delta_ts=[delta_t] + self.prev_delta_ts
		self.prev_time = t


def test_sep():
	n=Poly([5,1])
	d=Poly([5,11,1])
	p=Pade(n,d)
	s=separate(p)
	s.num.coeff=[r.real for r in s.num.coeff]
	s.denom.coeff=[r.real for r in s.denom.coeff]
	return s

sep=test_sep()
K=sep.num.coeff
s=sep.denom.coeff
def zeta(i,delta_n):
	return -s[i]*delta_n
def Phi(i,delta_n):
	return exp(s[i]*delta_n)
def q1(a:int,i:int,delta_n:float):
	return 0.
def q2(a:int,i:int,delta_n:float) -> float:
	phi = Phi(i,delta_n)
	if delta_n==0:
		return 0
	zi=zeta(i,delta_n)
	q0=(delta_n/zi**2)*(-1+zi+phi)
	q1=(delta_n/zi**2)*(1-(1+zi)*phi)
	return (q0,q1)[a]
def q3(a:int,i:int,delta_n:float) -> float:
	return 0.
def q4(a:int,i:int,delta_n:float) -> float:
	return 0.
def enqueue(val,lis):
	return [val] + lis[:-1]

def step(solv:Solver,t:float,v:float):
	solv.num_iter+=1
	x=solv.idxs.index('x')
	solv.update(t,v)
	dt=solv.prev_delta_ts
	final_val:float=0
	for i in range(len(K)):
		temp=Phi(i,dt[0])*solv.prev_values[i][0]
		for j in range(len(solv.K)):
			q=solv.q(j,i,dt[j])
			temp+=K[i]*q*solv.prev_values[x][0]
		final_val+=temp
		solv.prev_values[i]=enqueue(temp,solv.prev_values[i])
	y=solv.idxs.index('y')
	solv.prev_values[y]=enqueue(final_val,solv.prev_values[y])
	return final_val
	


if __name__=="__main__":
	K=[0.9478,0.0577]
	Kinv = [1/k for k in K]
	s=[-2.105,-0.095]
	sinv = [-_s for _s in s]
	solv=Solver(Kinv,sinv)
	i=0
	times=[]
	while i<40:
		i+=uniform(1e-5,1e-4)
		times.append(i)
	times=np.arange(0,40,1e-2)
	outputs=[]
	inputs=[]
	for i in times:
		v=np.cos(i)
		inputs.append(v)
		outputs.append(step(solv,i,v))
		# print(f"{i},{outputs[-1]}")
	# ex = np.loadtxt("sara_example.csv",dtype=float,delimiter=',')
	# c = np.loadtxt("results.csv",dtype=float,delimiter=',')
	# plt.plot(ex[:,0],ex[:,1])
	# plt.plot(c[:,0],c[:,1])
	plt.plot(times,np.cos(times))
	x_idx = solv.idxs.index('x')
	int_v = solv.prev_values[x_idx][::-1][5:]
	plt.plot(times,int_v)
	plt.show()
