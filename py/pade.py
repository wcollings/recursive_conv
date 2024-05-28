from typing import Sequence
from numpy import linspace, ones,zeros,chararray,array,logspace
from numpy.linalg import solve
from math import log10,ceil,exp,atan
from si_num import to_si
from matplotlib import pyplot as plt

def print_pade(aa:Sequence[float],bb:Sequence[float]):
	num=""

	# Create numerator string
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
	# num=straight_eval(aa,x)
	num=unroll(aa,x)
	# denom=straight_eval(bb,x)
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
	for i in range(N):
		for j in range(N):
			if i+M>=j:
				lower[i][j]="s"+str(i-j+M)
		y[i]="-s"+str(M+i+1)
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
	a=1e-9
	b=2.8e-9
	c=800e-9
	f0=2e4
	res=(0.6366*a)*atan(-c*(f-f0))+b;
	return res
def solve_system(s:list[float],M:int,N:int) -> tuple[tuple[float,...],tuple[float,...]]:
	"""
	Create a 1-point Pade approximation, given a list of coefficients and the length of the top and bottom polynomials
	"""
	y=zeros((N,1))
	lower=zeros((N,N))
	for i in range(N):
		for j in range(N):
			if i+M>=j:
				lower[i][j]=s[i-j+M]
		y[i]=-s[M+i+1]
	b=solve(lower,y)
	y=ones((M+1,1))
	upper=zeros((M+1,M+1))
	for i in range(M+1):
		for j in range(M+1):
			if i>=j:
				upper[i][j]=s[i-j]
	b_offset=ones((M+1,1))
	b_offset[1:]=b[:M]
	a=upper@b_offset
	a=(*map(float,a.flatten()),)
	b=(1,*map(float,b.flatten()),)
	return a,b
def get_err(aa,bb,plot=False):
	xmax=10000
	xs=logspace(1,8,10000)
	sig=tuple(map(L,xs))
	# sig=tuple(map(lambda x:1/(1+exp(x)),xs))
	sig=array(sig)
	approx=tuple(map(lambda x:eval_pade(aa,bb,x),xs))
	approx=array(approx)
	err=(sig-approx)/sig
	if plot:
		plt.semilogx(xs,sig)
		plt.semilogx(xs,approx)
		# plt.ylim(-2,2)
		plt.show()
	return sum(err)


def hf(n,N):
	num=3.14159-2*N*(-1)**n
	den=(2*n+1)
	return num/den

def lf(n):
	return 1/(2*n+1)
def twopoint(N):
	cs=[]
	for n in range(N):
		cs.append(hf(N+1-n,N))
		cs.append(0)
	for n in range(N):
		cs.append(lf(n))
		cs.append(0)
	cs.pop()
	b=[1]*N
	a=[1]*N
	for i in range(N):
		coeffs=[b[j]*cs[i-j+N] for j in range(N)]
		a[i]=sum(coeffs)
	print_pade(a,b)
	return a,b

if __name__=="__main__":
	maclaurin_coeffs=[2.158517e-09,-1.332268e-16,-2.467162e-19,1.827528e-21,-7.614698e-24,2.030586e-26,0.000000e+00,-2.558057e-31,]
	aa=[5.4284e34,6.9769e40,7.0833e46,3.8223e52][::-1]
	bb= [-4.8976e32, 2.708e43, 2.618e49, 2.731e55, 1.256e61][::-1]
	sample_matrix(3,4)
	aa,bb=solve_system(maclaurin_coeffs,3,4)
	print_pade(aa,bb)
	# aa,bb=twopoint(5)
	get_err(aa,bb,True)
	# for i in range(1,4):
	# 	for j in range(i+1,7-i+1):
	# 		aa,bb=solve_system(maclaurin_coeffs,i,j)
	# 		err=get_err(aa,bb,(i==1 and j==6))
	# 		print(f"M={i}\tN={j}\t{err=}")
	# plt.plot(xs,sig)
	# plt.plot(xs,approx)
	# plt.ylim(-2,2)
	# plt.show()
