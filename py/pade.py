from typing import Sequence
from ddiff import ddiff
from numpy import linspace, ones,zeros,chararray,array,logspace,ndarray
from numpy.linalg import solve
from math import log10,ceil,exp,atan
from si_num import to_si
from matplotlib import pyplot as plt

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
			output[i]+=cm[i+j][j]*aa[i+j]*(center**(j+1))
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
		plt.ylim((1.78e-9,3e-9))
		plt.show()
	return sum(err)


def hf(n,N):
	num=3.14159-2*N*(-1)**n
	den=(2*n+1)
	return num/den

def lf(n):
	return 1/(2*n+1)

if __name__=="__main__":
	start=int(1e6)
	step=20
	ne=8
	xs=list(range(start,start+(ne*step),step))
	ys=tuple(map(L,xs))
	dd=ddiff(xs,ys)
	# maclaurin_coeffs=[2.158517e-09,-1.332268e-16,-2.467162e-19,1.827528e-21,-7.614698e-24,2.030586e-26,0.000000e+00,-2.558057e-31,]
	# aa=[5.4284e34,6.9769e40,7.0833e46,3.8223e52][::-1]
	# bb= [-4.8976e32, 2.708e43, 2.618e49, 2.731e55, 1.256e61][::-1]
	# sample_matrix(3,4)
	aa,bb=solve_system(dd,3,4)
	aa=recenter_poly(aa,start)
	bb=recenter_poly(bb,start)
	print_pade(aa,bb,start)
	get_err(aa,bb,True)
