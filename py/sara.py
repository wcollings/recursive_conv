"""
Usage:
	pretty much have to hard-code in K,s, and the order of the state variables
	Call function `sara` at each time step, passing in the current sim time and the list of new state space variables it needs
	That in turn calls the specific solver method for each equation
"""
from math import exp, sin
from pade import Pade,separate,L,solve_system
from deriv import take_derivs
from poly import Poly
from matplotlib import pyplot as plt
from random import uniform
from typing import Callable
import numpy as np
from scipy.signal import find_peaks
import cmath
import pandas as pd

class Solver:
	prev_time:float
	prev_delta_ts:list[float]
	prev_values:list[list[complex]]
	raw_integral:list[complex]
	prev_x:float
	num_iter:int
	tstep:float
	K:list
	s:list
	idxs:list[str]
	k0:float
	q:Callable
	def __init__(self,p:Pade):
		self.prev_x=0
		self.tstep=-1
		self.k0=p.k0
		self.prev_time=0
		self.num_iter=0
		self.K=p.num.coeff
		self.s=p.denom.coeff
		self.prev_delta_ts=[0]*8
		self.raw_integral=[0]*3
		self.prev_values=[[0]*5 for i in range(len(p.num.coeff)+2)]
		self.idxs=['y'+str(i) for i in range(len(p.num.coeff))] + ['y','x']
		if len(p.num.coeff)==1:
			self.q=q1
		elif len(p.num.coeff)==2:
			self.q=q2
		elif len(p.num.coeff)==3:
			self.q=q3
		else:
			self.q=q4

	def __repr__(self):
		s=f"""{self.num_iter} iterations. Equation was:
			K0={self.k0}
			K={self.K}
			s={self.s}
			"""
		return s
	def update(self,t,x):
		x_idx = self.idxs.index('x')
		prev_x = self.prev_values[x_idx][0]
		delta_t=self.tstep
		if self.tstep==-1:
			delta_t = t-self.prev_time
		new_x = prev_x + (delta_t*(x+self.raw_integral[0])/2)
		self.raw_integral.pop()
		self.raw_integral = [new_x] + self.raw_integral
		update_pat = [4,2,1]
		new_integral = sum(map(lambda x:x[0]*x[1],zip(update_pat,self.raw_integral)))/7
		# new_integral = sum(self.raw_integral)/len(self.raw_integral)
		self.prev_values[x_idx].pop()
		self.prev_values[x_idx]= [new_integral] + self.prev_values[x_idx]
		if len(self.prev_values[x_idx]) > 4:
			print("List is expanding!")
			self.prev_values[x_idx]=self.prev_values[x_idx][:4]
		self.prev_delta_ts=[delta_t] + self.prev_delta_ts
		self.prev_time = t
		self.prev_x=x
		return new_x


def test_sep():
	n=Poly([5,1])
	d=Poly([5,11,1])
	p=Pade(n,d)
	s=separate(p)
	s.num.coeff=[r.real for r in s.num.coeff]
	s.denom.coeff=[r.real for r in s.denom.coeff]
	return s

# sep=test_sep()
# K=sep.num.coeff
# s=sep.denom.coeff
def zeta(si,delta_n):
	return -si*delta_n
def Phi(si,delta_n):
	return cmath.exp(si*delta_n)
def q1(a:int,i:int,delta_n:float):
	return delta_n
def q2(a:int,i:float,delta_n:float) -> complex:
	# if delta_n==0:
	# 	return 0
	# zi=zeta(i,delta_n)
	# phi = Phi(i,delta_n)
	# q0=(delta_n/zi**2)*(-1+zi+phi)
	# q1=(delta_n/zi**2)*(1-(1+zi)*phi)
	q0=0
	q1=delta_n
	return (q0,q1)[a]
def q3(a:int,i:int,delta_n:float) -> float:
	zi=zeta(i,delta_n)
	phi=Phi(i,delta_n)
	# q0=delta_n/(2*zi**3)*(2-3*zi+2
def q4(a:int,i:int,delta_n:float) -> float:
	return 0.
def enqueue(val,lis):
	return [val] + lis[:-1]

def step(solv:Solver,t:float,v:float):
	solv.num_iter+=1
	x=solv.idxs.index('x')
	y=solv.idxs.index('y')
	curr_int=solv.update(t,v)
	dt=solv.prev_delta_ts
	# if abs(dt[1]-dt[0])/dt[1] > 0.75:
	# 	print("Sudden change in delta_t")
	# 	print(f"{t=}")
	final_val:complex=solv.k0*curr_int
	num_terms = len(solv.K)
	for i in range(num_terms):
		temp=Phi(solv.s[i],dt[0])*solv.prev_values[i][0]
		for j in range(len(solv.K)):
			q=solv.q(j,solv.s[i],dt[j])
			temp+=solv.K[i]*q*solv.prev_values[x][j]
		final_val+=temp
		solv.prev_values[i]=enqueue(temp,solv.prev_values[i])
	solv.prev_values[y]=enqueue(final_val,solv.prev_values[y])
	return final_val.real
	

def load_approx(i)->Pade:
	if i==0:
		r1=1
		r2=1
		c1=1e-6
		c2=5e-6
		p=separate(Pade(Poly([r2*c2,1]),Poly([r1*c1*r2*c2,r1*c1+r1*c2+r2*c2,1])))
		return p
	if i==1:
		x0=1e6
		x1=1e9
		l=lambda s:L(s)-L(x1)
		derivs=take_derivs(l,x0,3,40)
		# print(derivs)
		p1=Poly(derivs[::-1]).recenter(x0)
		p=solve_system(p1,1,2)
		# p.k0=1e8
		p.k0=L(x1)
		print(f"K0={L(x1)} -> {1/L(x1)}H")
		return separate(p)
	return Pade(Poly([0]),Poly([0]))

def run(solv, fname:str):
	datafile=pd.read_csv(fname)
	datafile['i']=0
	datafile['int_v']=0
	datafile['err']=0
	x=solv.idxs.index('x')
	for i in range(datafile.shape[0]):
		v=datafile['v'][i]
		t=datafile['t'][i]
		res=step(solv,t,v)
		datafile['i'][i]=res
		known=datafile['known_i'][i]
		datafile['err'][i]=abs(known-res)/res
		datafile['int_v'][i]=solv.prev_values[x][0]
	return datafile



if __name__=="__main__":
	p=load_approx(1)
	print(p)
	# K=[k.real for k in p.num.coeff]
	# s=[sig.real for sig in p.denom.coeff]
	# with open("../tcl/coeff.csv","w") as fp:
	# 	fp.write(p.get_coeff_str())
	solv=Solver(p)
	pi=3.14
	end=1000
	freq=1e7
	tstep=20/freq

	solv=Solver(p)
	# solv.k0*=1.4
	# solv.K = [K*2 for K in solv.K]
	f_1e5=run(solv,"vsig_low.csv")
	solv=Solver(p)
	f_1e8=run(solv,"vsig.csv")
	fig,axs = plt.subplots(2,2,sharex='row') #pyright:ignore
	axs[0,0].plot(f_1e5.t,f_1e5.i,label="predicted")
	axs[0,0].plot(f_1e5.t,f_1e5.known_i, label="static L, 1.3nH")
	axs[0,0].grid(True)
	axs[0,0].set_title("f=100kHz")
	axs[0,0].legend()
	axs[0,1].plot(f_1e5.t,f_1e5.err)
	axs[0,1].grid(True)
	axs[0,1].set_title("Error")

	axs[1,0].set_title("f=10MHz")
	axs[1,0].plot(f_1e8.t,f_1e8.i,label="predicted")
	axs[1,0].plot(f_1e8.t,f_1e8.known_i,label="static L, 1.3nH")
	axs[1,0].grid(True)
	axs[1,0].legend()
	axs[1,1].plot(f_1e8.t,f_1e8.err)
	axs[1,1].grid(True)
	axs[1,1].set_title("Error")
	# axs[0].plot(f_1e6.t,f_1e6.int_v)
	# axs[0].grid(True)
	plt.show()
