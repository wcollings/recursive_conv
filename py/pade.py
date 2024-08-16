from typing import Sequence, Iterable
from ddiff import ddiff
from poly import Poly,synth_div
from numpy import linspace, ones,zeros,chararray,array,logspace,ndarray,vectorize
import numpy as np
from numpy.linalg import solve
from math import log10,ceil,exp,atan, factorial as fac
from deriv import take_derivs
from matplotlib import pyplot as plt

class Pade:
	num:Poly
	denom:Poly
	sep:bool
	k0:float
	def __init__(self,n:Poly,d:Poly,k0:float=0.):
		self.num=n
		self.denom=d
		self.sep=False
		self.k0=k0
	def __call__(self,x:float) -> float:
		if self.sep:
			ret=self.k0
			for K,s in zip(self.num.coeff,self.denom.coeff):
				ret += K/(x-s)
			return ret
		return self.k0 + (self.num(x)/self.denom(x))
	def __repr__(self):
		out=''
		if self.sep:
			center="r(x)= " + str(self.k0) + " + "
			num = " "*len(center)
			denom = " "*len(center)
			temp_n =""
			temp_d=""
			for n,d in zip(self.num.coeff,self.denom.coeff):
				temp_n = f"{n:1.2e}"
				temp_d = f"(x{d:-1.2e})"
				if len(temp_n) > len(temp_d):
					center+="+ " + "-"*len(temp_n)
					num+="  " + temp_n
					denom +="  " + f"{temp_d: ^{len(temp_n)}}"
				else:
					center+=" + " + "-"*len(temp_d)
					num+="   " + f"{temp_n: ^{len(temp_d)}}"
					denom +="   " + temp_d
			out = num+"\n" + center+"\n" +denom
		else:
			num=str(self.num)
			denom=str(self.denom)
			str_len=max(map(len,[num,denom]))
			center_len = len("r(x)= " + str(self.k0) + " + ")
			center="r(x)= " + str(self.k0) + " + " + "-"*len(max(num,denom))+'\n'
			out=f"{' '*center_len}{num:^{str_len}s}\n"+center+f"{' '*center_len}{denom:^{str_len}s}"
		return out

def pascal(row,invert=True) -> ndarray:
	if row==1:
		return ones(1)
	res = zeros((row,row))
	res[:,0]=1
	for j in range(1,row):
		for k in range(j,0,-1):
			res[j,k]=res[j-1,k]+res[j-1,k-1]
	if invert:
		for i in range(row):
			for j in range(row-1):
				res[i,j] *= (-1)**(j)
	return res
	
def recenter_poly(aa:Sequence[float],center:float):
	cm = pascal(len(aa)) # Coefiencient Multipliers 
	output=[0.]*len(aa)
	for i in range(len(aa)):
		for j in range(len(aa)-i):
			M=cm[i+j][j]
			an=aa[i+j]
			output[i]+=cm[i+j][j]*aa[i+j]*(center**(j))
	return output

def print_pade(aa:Sequence[float],bb:Sequence[float], base:float=0):
	num=""
	# Create numerator string
	if base != 0:
		aa=recenter_poly(aa,base)
		bb=recenter_poly(bb,base)
	for i,_a in enumerate(aa):
		a=to_si(_a,sign=True)
		if i==0:
			a=to_si(_a,sign=False)
			num+=str(a)
		elif i==1:
			num+=str(a) + "x"
		else:
			num+=str(a) + "x^" + str(i)

	denom=""
	bb=(1.,*bb)
	for i,_b in enumerate(bb):
		b=to_si(_b,sign=True)
		pow="x^"+str(i)
		if _b==1:
			b="+"
		elif _b==-1:
			b="-"
		if i==0:
			b=to_si(_b,sign=False)
			pow=""
		elif i==1:
			pow="x"
		denom+=b+pow
	str_len=max(map(len,[num,denom]))
	center="r(x)= "+ "-"*(str_len)
	print(f"{' '*6}{num:^{str_len}s}")
	print(center)
	print(f"{' '*6}{denom:^{str_len}s}")

def eval_pade(aa,bb, x):
	num=unroll(aa,x)
	denom=unroll(bb,x)
	return num/denom

def unroll(aa:list[float],x:float):
	if len(aa)==1:
		return aa[0]
	return aa[0]+(unroll(aa[1:],x)*x)

def straight_eval(aa:list[float],x:float) -> float:
	res=0
	for i,v in enumerate(aa):
		res+=(v*x**(i))
	return res

def sample_matrix(M,N):
	"""
	prints the coefficients used for solving the one-point pade approximation in index
	form. Useful for debugging, but won't do much besides that
	"""
	# Solve for Bs
	charlen=ceil(log10(M+N+1))+1
	y=chararray((N,1),itemsize=charlen+1)
	lower=chararray((N,N),itemsize=charlen)
	lower[:]='0'*charlen
	b_arr = ["s"+str(i) for i in range(M+N+1)]
	for i in range(N):
		for j in range(N):
			if i+M>=j:
				lower[i][j]=b_arr[i-j+M]
				# lower[i][j]="s"+str(i-j+M)
		y[i]="-"+b_arr[M+i+1]
		# y[i]="-s"+str(M+i+1)
	for r in range(N):
		print('['+','.join(lower[r].decode()),end=']')
		print(y[r].decode())
	print()
	# Solve for As
	y=chararray((M+1,1),itemsize=charlen+1)
	b=chararray((M+N+1,1),itemsize=charlen)
	b[:]='0'*(charlen-1)+'1'
	A=chararray((M+1,M+1),itemsize=charlen)
	A[:]='0'*charlen
	for i in range(M+1):
		y[i]="a"+str(i)
	for i in range(1,M+1):
		b[i]="b"+str(i)
	# for i in range(M+1,M+N+1):
	for i in range(M+1):
		for j in range(M+1):
			if i>=j:
				A[i][j]="s"+str(i-j)
	for r in range(M+1):
		print('['+','.join(A[r].decode()),end=']')
		print(b[r].decode(),end='=')
		print(y[r].decode())

def L(f):
	a=21.81e-9
	b=25.35e-9
	c=4.01e-7
	f0=374.06e3
	res=(0.6366*a)*atan(-c*(f-f0))+b;
	return 1/res
def solve_system(p:Poly,M:int,N:int) -> Pade:
	"""
	Create a Pade approximation, given a Taylor series and the length of the top and bottom polynomials
	len(p) must equal (M+N+1)
	"""
	print(p.coeff)
	s=p.coeff[::-1]
	# print("s")
	# print(s)
	y=zeros((N,1))
	lower=zeros((N,N))
	for i in range(N):
		for j in range(N):
			if i+M>=j:
				lower[i][j]=s[i-j+M]
		print(f"s[{M}+{i}+1]={s[M+i+1]}")
		y[i]=-s[M+i+1]
	# print(lower)
	# print(y)
	b=solve(lower,y)
	# print(b)
	y=ones((M+1,1))
	upper=zeros((M+1,M+1))
	for i in range(M+1):
		for j in range(M+1):
			if i>=j:
				upper[i][j]=s[i-j]
	b_offset=ones((M+1,1))
	b_offset[1:]=b[:M]
	# print()
	# print(upper)
	a=upper@b_offset
	a=list(map(float,a.flatten()))
	b=[1]+list(map(float,b.flatten()))
	return Pade(Poly(a[::-1]),Poly(b[::-1]))
	# return Pade(Poly(a[::-1]).recenter(-2e6),Poly(b[::-1]).recenter(-2e6))

def get_err(rep,plot=False):
	xmax=10000
	xs=logspace(1,9,10000)
	sig=tuple(map(L,xs))
	if not isinstance(rep,Iterable):
		rep=[rep]
	# sig=tuple(map(lambda x:1/(1+exp(x)),xs))
	sig=array(sig)
	approx=[]
	for r in rep:
		# print(r.disc)
		approx.append(vectorize(r)(xs))
		# err=(sig-approx[-1])/sig
	if plot:
		plt.semilogx(xs,sig,label="actual")
		for n,i in enumerate(approx):
			plt.semilogx(xs,i,label=f"approx #{n}")
		# plt.ylim((1.4e-9,4.55e-9))
		plt.legend()
		plt.show()
	# return sum(err)

def separate(rep:Pade):
	mat = zeros((len(rep.denom)-1,len(rep.denom)-1),dtype=np.complex64)
	# print(mat)
	roots = rep.denom.get_roots()
	# print(f"The roots of {rep.denom} are:\n {roots}")
	K=[]
	for r in roots.coeff:
		temp = synth_div(rep.denom,r)
		K.append(rep.num(r)/temp(r))
	# print(K)
	out=Pade(Poly(K),roots) #pyright:ignore
	out.k0=rep.k0
	out.sep=True
	return out

def twopt(rep,start):
	def inner(x):
		return L(start) + rep(x)
	inner.disc = f"{start} + ({str(rep)})"
	return inner

if __name__=="__main__":
	x0=1e6
	x1=1e9
	l=lambda s:L(s)-L(x1)
	derivs=take_derivs(l,x0,3,40)
	# p1=Poly(derivs[::-1]).recenter(x0)
	# derivs=list(map(float,open("L_deriv.csv").readline().split(",")))[:4]
	poly=Poly(derivs[::-1]).recenter(x0)
	# print(poly)
	p=solve_system(poly,1,2)
	p.k0=L(x1)
	s=separate(p)
	print(s)
	# get_err([p,s],True)
