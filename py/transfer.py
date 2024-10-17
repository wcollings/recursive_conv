from typing import Callable
from pade import Pade
from poly import Poly
import numpy as np
import cmath

class Transfer:
	Ki:list[complex]
	si:list[complex]
	K0:complex
	Kn:complex
	def __init__(self,K0,Ki,si,d) -> None:
		self.K0=K0
		self.Ki=Ki
		self.si=si
		self.Kn=d

def zeta(si,delta_n):
	return -si*delta_n
def Phi(si,delta_n):
	return cmath.exp(si*delta_n)
def q1(a:int,i:int,delta_n:float):
	return delta_n
def q2(a:int,i:float,delta_n:float) -> complex:
	zi=zeta(i,delta_n)
	q0=delta_n/2
	q1=(delta_n/2)*(1-zi)
	return (q0,q1)[a]
def q3(a:int,i:int,delta_n:float) -> float:
	zi=zeta(i,delta_n)
	phi=Phi(i,delta_n)
	return 0.
def q4(a:int,i:int,delta_n:float) -> float:
	return 0.

def enqueue(val,lis):
	return [val] + lis[:-1]

class TSolver:
	prev_time:float
	tn:list[float]
	xn:list[float]
	dxn:list[complex]
	yn:list[complex]
	eq:Transfer
	order:int
	num_iter:int
	q:Callable
	def __init__(self,order,eq,ic=0) -> None:
		self.prev_time=0
		self.tn=[0]*order
		self.xn=[0]*order
		self.dxn=[0]*order
		self.yn=[0]*order
		self.eq=eq
		self.order=order
		self.num_iter=0
		self.q=[q1,q2,q3,q4][order]
		if (ic !=0):
			self.tn=[1e-9]*order
			self.xn=[ic]*order
			last_res=0
			while (abs(res:=self.step(1e-9,ic)-last_res) > 1e-5):
				self.prev_time=-1e-9
				last_res=res

	def step(self,t:float,x:float):
		self.num_iter+=1
		dxn = x + self.dxn[0] - self.xn[0]
		self.xn=enqueue(x,self.xn)
		self.tn=enqueue(t-self.prev_time,self.tn)
		self.dxn=enqueue(dxn,self.dxn)
		final_val:complex=self.eq.K0*x + self.eq.Kn*dxn
		num_terms = len(self.eq.Ki)
		for i in range(num_terms):
			temp=Phi(self.eq.si[i],self.tn[0])*self.yn[i]
			for j in range(len(self.eq.Ki)):
				q=self.q(j,self.eq.si[i],self.tn[j])
				temp+=self.eq.Ki[i]*q*self.xn[j]
			final_val+=temp
			self.yn[i]=temp
		self.prev_time=t
		return final_val.real
