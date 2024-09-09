from si_num import from_si,to_si
from typing import Sequence, Iterable
from poly import Poly,synth_div
from numpy import ones,zeros,chararray,array,logspace,ndarray,vectorize
import numpy as np
from numpy.linalg import solve
from math import log10,ceil,exp,atan, factorial as fac
from deriv import take_derivs
from matplotlib import pyplot as plt
from plotting import figure_wrapper
import pandas as pd
import cmath

import warnings
warnings.filterwarnings("ignore")
class max_monad:
	val:complex
	highest:float
	def __init__(self, val):
		self.val=val
		if isinstance(val,complex):
			val = max(val.real,val.imag)
		self.highest=val
	def bind(self,val) -> 'max_monad':
		orig=val
		if isinstance(val,complex):
			val = max(val.real,val.imag)
		if val > self.highest:
			return max_monad(orig)
		return self
	def __repr__(self):
		return f"Highest power was: {int(log10(self.highest))} ({self.val:1.2e})"

def complex_str(a:complex|float):
	ic=lambda a: isinstance(a,complex) and a.imag!=0
	if ic(a):
		return f"{a.real:-f} {a.imag:-f} "
	return f"{a:-f} 0 "

def format_float(f:complex|float) -> str:
	if isinstance(f,float):
		rep=f"{f:+1.3e}"
		man,exp=rep.split("e")
		if exp.startswith("+"):
			exp=exp[1:]
		exp=exp.strip("0")
		return man+r"\text{e}"+exp
	real=format_float(f.real)
	if f.imag>0:
		sign="+"
	else:
		sign="-"
	imag=sign+"j"+format_float(f.imag)[1:]
	return real+imag


class Pade:
	num:Poly
	denom:Poly
	sep:bool
	freq_domain:bool
	k0:float
	def __init__(self,n:Poly,d:Poly,k0:float=0.):
		self.freq_domain=True
		self.num=n
		self.denom=d
		self.sep=False
		self.k0=k0
	def __call__(self,x:float) -> float:
		if not self.freq_domain:
			ret=self.k0*x
			for K,s in zip(self.num.coeff,self.denom.coeff):
				this_term = K*cmath.exp(s*x)
				ret+=this_term
			return ret
		if self.sep:
			ret=self.k0
			for K,s in zip(self.num.coeff,self.denom.coeff):
				this_term=K/(x-s)
				ret += this_term
			return ret
		elif not self.freq_domain:
			ret=self.k0*x
			for K,s in zip(self.num.coeff,self.denom.coeff):
				this_term = K*exp(s*x)
				ret+=this_term
			return ret
		return self.k0 + (self.num(x)/self.denom(x))
	def __repr__(self):
		out=''
		if not self.freq_domain:
			out=f"r(x) = {self.k0} "
			for k,s in zip(self.num.coeff,self.denom.coeff):
				out += f" + ({k})exp({s}t)"
			return out
		if self.sep:
			center=f"r(x) = {self.k0:e} " 
			num = " "*len(center)
			denom = " "*len(center)
			temp_n =""
			temp_d=""
			for n,d in zip(self.num.coeff,self.denom.coeff):
				temp_n = f"{n:1.2e}"
				temp_d = f"(x{-d:+1.2e})"
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
			center=f"r(x) = {self.k0:e} + " 
			center_len = len(center)
			center+="-"*len(max(num,denom))+'\n'
			out=f"{' '*center_len}{num:^{str_len}s}\n"+center+f"{' '*center_len}{denom:^{str_len}s}"
		return out
	def get_coeff_str(self)->str:
		s=complex_str(self.k0)
		for coeff in self.num.coeff:
			s+=complex_str(coeff)
		if len(self.num.coeff) < 4:
			s+= "0 0 "*(4-len(self.num.coeff))
		for coeff in self.denom.coeff:
			s+=complex_str(coeff)
		if len(self.denom.coeff) < 4:
			s+= "0 0 "*(4-len(self.denom.coeff))
		return s
	def highest_power(self):
		found=max_monad(self.k0)
		for k in self.num.coeff:
			found=found.bind(k)
		return found
		# print(f"The highest power I found was {int(log10(found.highest))} ({found.val:1.3e})")
	def to_latex(self):
		s="H(s)="
		if not self.freq_domain:
			s="h(t)="
			s+=format_float(self.k0)
			for n,d in zip(self.num.coeff,self.denom.coeff):
				s+=r"+("+format_float(n)+r")\exp(("+format_float(-d)+")t)"
			return s
		if self.sep:
			if self.k0:
				s+=format_float(self.k0)
			for n,d in zip(self.num.coeff,self.denom.coeff):
				s+=r"+\frac{"+format_float(n)+r"}{s"+format_float(-d)+"}"
			return s
		return ""


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
	# print(p.coeff)
	s=p.coeff[::-1]
	# print("s")
	# print(s)
	y=zeros((N,1))
	lower=zeros((N,N))
	for i in range(N):
		for j in range(N):
			if i+M>=j:
				lower[i][j]=s[i-j+M]
		# print(f"s[{M}+{i}+1]={s[M+i+1]}")
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
	xs=logspace(1,9,10000)
	sig=tuple(map(L,xs))
	if not isinstance(rep,Iterable):
		rep=[rep]
	sig=array(sig)
	approx=[]
	err=[]
	saber_overlay=pd.read_csv("lf_matching.csv")
	for r in rep:
		r.freq_domain=True
		res = vectorize(r)(xs)
		approx.append(vectorize(r)(xs))
		err.append(np.sum(np.divide(np.abs(res-sig),sig)))
	if plot:
		with figure_wrapper(outf="L_compare.png",show=True) as fw:
			fw.plot2(xs,1/sig,sig,xlab="Frequency (Hz)",ylab1="L(f)",ylab2="1/L(f)",name="Original function")
			linv = [1/l for l in saber_overlay.l]
			fw.plot2(saber_overlay.f,saber_overlay.l,linv, name="From Saber")
			fw.plot2(xs,1/approx[0],approx[0],xlab="Frequency (Hz)",ylab1="L(f)",ylab2="$L^{-1}(f)$",name="Approximation")
			fw.set_fontsize(15)
	return err

def separate(rep:Pade):
	# mat = zeros((len(rep.denom)-1,len(rep.denom)-1),dtype=np.complex64)
	roots = rep.denom.get_roots()
	K=[]
	for r in roots.coeff:
		temp = synth_div(rep.denom,r)
		K.append(rep.num(r)/temp(r))
	out=Pade(Poly(K),roots) #pyright:ignore
	out.k0=rep.k0
	out.sep=True
	return out

def twopt(rep,start):
	def inner(x):
		return L(start) + rep(x)
	inner.disc = f"{start} + ({str(rep)})"
	return inner
def create_approx(x):
		x0=x*1e6
		x1=1e9
		l=lambda s:L(s)-L(x1)
		derivs=take_derivs(l,x0,3,40)
		p1=Poly(derivs[::-1]).recenter(x0)
		p=solve_system(p1,1,2)
		# print(p)
		s=separate(p)
		# print(s)
		s.k0=L(x1)
		s.sep=True
		return s,p1

def time_domain(rep:Pade):
	rep.freq_domain=False
	rep.num.coeff = [n*100 for n in rep.num.coeff]
	# rep.denom.coeff = [n*100 for n in rep.denom.coeff]
	print(rep)
	with open("l_t_latex.txt",'w') as f:
		f.write(rep.to_latex())
		f.write("\n")
		rep.freq_domain = True
		f.write(rep.to_latex())
		rep.freq_domain=False
	ts=np.arange(0,1e-4,1e-8)
	xs=1/np.vectorize(rep)(ts)
	ts=ts*1e6
	xs=xs*1e3
	with figure_wrapper("L_time_domain.png", show=True) as fw:
		fw.plot(ts,xs)
		fw.set_xlim(0,1e2)
		fw.set_ylim(-6,2.6)
		fw.set_labels("Time (us)","Inductance (mH)")
		fw.set_title("Time-Domain Inductance response (100K)")
		fw.set_fontsize(15)

if __name__=="__main__":
	min_err=1e9
	min_x0=1e9
	approxs=[]
	best_approx=Pade(Poly([]),Poly([]))
	xs = np.linspace(1,12,12)
	df=pd.DataFrame(0,index=xs,columns=['err','pow'])
	s,p=create_approx(1)
	best_approx=s
	# print(best_approx.to_latex())
	# s.num.coeff = [k*2 for k in s.num.coeff]
	test_f = [1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10]
	fbase=np.linspace(1,10,10)
	f1=fbase*1e5
	f2=fbase*1e6
	fs=f1.tolist() + f2.tolist()
	ls=tuple(map(lambda s:1/best_approx(s),fs))
	ls = [1e9/best_approx(f).real for f in fs]
	pd.DataFrame(data={'f':fs,'l':ls}).to_csv("real_lf.csv")
	df=pd.DataFrame(0,index=test_f,columns=['L','H','delta'])
	for freq in test_f:
		a=1/L(freq)
		b=1/best_approx(freq).real
		df.at[freq,'L']=a
		df.at[freq,'H']=b
		df.at[freq,'delta']=abs(a-b)/a*100
	# print(df)
	time_domain(best_approx)
	# get_err([best_approx],True)
