from math import sin,cos,pi
from scipy.integrate import cumulative_trapezoid as cumtrapz
from random import uniform
import numpy as np
from matplotlib import pyplot as plt
from poly import Poly
from numpy.linalg import solve
import numpy as np

def cubic_spline(ts,xs,i):
	t2=ts[i]
	t1=ts[i-1]
	t0=ts[i-2]
	x2=xs[i]
	x1=xs[i-1]
	x0=xs[i-2]
	delta_1=(t1-t0)
	delta_2=(t2-t1)
	mat=np.zeros((5,5))
	b=np.zeros((5,1))
	mat[0,0]=delta_1
	mat[0,1]=delta_1**2
	b[0]=x1-x0

	mat[1,0]=1
	mat[1,1]=3*delta_1**2
	mat[1,2]=-1

	mat[2,1]=6*delta_1
	mat[2,3]=-2
	
	mat[3,2]=delta_2
	mat[3,3]=delta_2**2
	mat[3,4]=delta_2**3
	b[3]=x2-x1

	mat[4,3]=2
	mat[4,4]=6*delta_2

	out=solve(mat,b)
	a0=x0
	b0=out[0]
	c0=0
	d0=out[1]

	a1=x1
	b1=out[2]
	c1=out[3]
	d1=out[4]
	s0=Poly([d0,c0,b0,a0]).recenter(t0)
	s1=Poly([d1,c1,b1,a1]).recenter(t1)
	return s0,s1


def progbar(prog,total):
	perc=int(100*(prog/float(total)))
	pperc=int(25*(prog/float(total)))
	bar = '=' *(pperc-1)+'>' + ' '*(25-pperc)
	print(f"\r[{bar}] {perc:2.1f}%",end='\r')

def L0(ts,xs,i):
	t0=ts[i-2]
	t1=ts[i-1]
	t2=ts[i]
	x0=xs[i-2]
	b=t2*t1
	s=(t2+t1)/2
	d=(t2-t1)/2
	c =x0/(((t0-s)**2)*(1-(d/(t0-s))**2))
	term1=(t2**3 - t1**3)/3
	term2=s*(t2**2-t1**2)
	term3=b*(t2-t1)
	ret= c*(term1-term2+term3)
	return ret

def L1(ts,xs,i):
	t0=ts[i-2]
	t1=ts[i-1]
	t2=ts[i]
	x1=xs[i-1]
	b=t2*t0
	s=(t2+t0)/2
	d=(t2-t0)/2
	c = x1/(((t1-s)**2)*(1-(d/(t1-s))**2))
	term1=(t2**3 - t1**3)/3
	term2=s*(t2**2-t1**2)
	term3=b*(t2-t1)
	ret= c*(term1-term2+term3)
	return ret

def L2(ts,xs,i):
	t0=ts[i-2]
	t1=ts[i-1]
	t2=ts[i]
	x2=xs[i]
	b=t1*t0
	s=(t1+t0)/2
	d=(t1-t0)/2
	c = x2/(((t2-s)**2)*(1-(d/(t2-s))**2))
	term1=(t2**3 - t1**3)/3
	term2=s*(t2**2-t1**2)
	term3=b*(t2-t1)
	ret= c*(term1-term2+term3)
	return ret

def lagrange_quad(ts,xs,i):
	return L0(ts,xs,i) + L1(ts,xs,i) + L2(ts,xs,i)
def lagrange_approx(ts,xs,i):
	t2=ts[i]
	t1=ts[i-1]
	t0=ts[i-2]
	x2=xs[i]
	x1=xs[i-1]
	x0=xs[i-2]
	def inner(t):
		l0 = (t-t1)*(t-t2)*x0/(t0-t1)/(t0-t2)
		l1 = (t-t0)*(t-t2)*x1/(t1-t0)/(t1-t2)
		l2 = (t-t0)*(t-t1)*x2/(t2-t0)/(t2-t1)
		return l0+l1+l2
	return inner

def trap(ts,xs,i):
	return (ts[i]-ts[i-1])*(xs[i]+xs[i-1])/2

def plot_lagrange(ts,xs):
	plt.plot(ts,xs)
	num_iter=len(xs)
	plt.plot(ts,xs)
	# plt.ion()
	for i in range(3,num_iter):
		t=ts[i-3:i]
		x=xs[i-3:i]
		# print(f"Interpolating over {t}")
		t_plot=np.linspace(ts[i-3],ts[i],10)
		s0,s1=cubic_spline(t,x,2)
		plt.plot(t_plot,s1(t_plot))
		# breakpoint()
	plt.show()

def create_input_arrays(num_pts):
	ts=np.cumsum(np.random.rand(num_pts,1)/1e9)
	xs=np.zeros(ts.shape)
	ys=np.zeros(ts.shape)
	t0=5e-8
	idx0=np.argmin(np.abs(ts-t0))
	idx1=np.argmin(np.abs(ts-2e-5))
	xs[idx0]=1
	ys[idx0]=xs[idx0]*(ts[idx0]-ts[idx0-1])
	xs[idx0+1]=2
	ys[idx0+1]=xs[idx0+1]*(ts[idx0+1]-ts[idx0])+ys[idx0]
	m=1e5
	t1=ts[idx0+2]
	for i in range(idx0+2,idx1):
		xs[i]=(ts[i]-t1)*m
		ys[i]=m*(ts[i]**2-t1**2)/2 #+ys[idx0+1]
	t2=ts[idx1]
	a=-1e5
	b=1e7
	for i in range(idx1,ts.shape[0]):
		xs[i]=0.1*np.cos(1e7*(ts[i]-t2))*np.exp((ts[i]-t2)*a)
		exp_term=np.exp((ts[i]-t2)*a)
		cos_term=a*np.cos(b*(ts[i]-t2))
		sin_term=b*np.sin(b*(ts[i]-t2))
		denom=a**2+b**2
		ys[i]=0.1*exp_term*(cos_term+sin_term)/denom + ys[idx1-1]

	return ts,xs,ys
def run_integral_compare(numpts:int,plot=False) -> tuple[int,int]:
	# ts=np.zeros((numpts,1))
	# xs=np.zeros((numpts,1))
	ts,xs,ys=create_input_arrays(numpts)
	# ys=np.zeros((numpts,1))
	yL=np.zeros((numpts,1))
	yC=np.zeros((numpts,1))
	yT=np.zeros((numpts,1))
	eL=np.zeros((numpts,1))
	eC=np.zeros((numpts,1))
	eT=np.zeros((numpts,1))
	teL=0
	teT=0
	teC=0
	# ts[1]=ts[0]+uniform(1e-10,1e-8)
	# xs[1]=np.cos(ts[1]*1e7)
	# ys[1]=np.sin(ts[1]*1e7)/1e7
	# ts[2]=ts[1]+uniform(1e-10,1e-8)
	# xs[2]=np.cos(ts[2]*1e7)
	# ys[2]=np.sin(ts[2]*1e7)/1e7
	ys[0]=0
	ys[1]=0
	for i in range(2,numpts):
		# ts[i]=ts[i-1]+uniform(1e-10,1e-8)
		# xs[i]=np.cos(ts[i]*1e7)
		# ys[i]=np.sin(ts[i]*1e7)/1e7
		# if (i>250 and i < 255):
		# 	xs[i]+=10
		# 	ys[i]+=10*(ts[i]-ts[250])
		# if i>= 255:
		# ys[i]+=10*(ts[255]-ts[250])
		yL[i]=lagrange_quad(ts,xs,i)+yL[i-1]
		_,s1=cubic_spline(ts,xs,i)
		s1=s1.integ()
		yC[i]=s1(ts[i])-s1(ts[i-1])+yC[i-1]
		yT[i]=trap(ts,xs,i)+yT[i-1]
		eL[i]=abs(yL[i]-ys[i])
		eC[i]=abs(yC[i]-ys[i])
		eT[i]=abs(yT[i]-ys[i])
		teL+=eL[i]
		teC+=eC[i]
		teT+=eT[i]
	if plot:
		fig,axs=plt.subplots(3,1,sharex=True)
		axs[1].plot(ts,yL,label="Lagrange integral")
		axs[1].plot(ts,yC,label="Cubic Spline integral")
		axs[1].plot(ts,yT,label="Trapezoidal method")
		axs[1].plot(ts,ys,label="Actual soln")
		axs[1].set_ylabel("Integral of signal")
		axs[1].legend()
		axs[2].plot(ts,eL,label="Lagrange error")
		axs[2].plot(ts,eC,label="Cubic Spline error")
		axs[2].plot(ts,eT,label="Trapezoidal error")
		axs[2].set_ylabel("Error of integrals")
		axs[2].legend()
		axs[0].plot(ts,xs)
		axs[0].set_ylabel("Signal")
		print(f"total Lagrange error: {teL}")
		print(f"total Spline error: {teC}")
		print(f"total Trapezoid error: {teT}")
		plt.show()
	return teL,teT

		

if __name__=="__main__":
	# ts=np.cumsum(np.random.rand(20,1)/1e-9)
	# xs=np.cos(ts*1e7)
	# plot_lagrange(ts,xs)
	numpts=100000
	num_runs=100
	err_L=np.zeros((num_runs,1))
	err_T=np.zeros((num_runs,1))
	# for i in range(num_runs):
	# 	progbar(i,num_runs)
	# 	err_L[i],err_T[i]=run_integral_compare(numpts)
	# print(f"Average error of Lagrange: {np.sum(err_L)/100}")
	# print(f"Average error of Trapezoid: {np.sum(err_T)/100}")
	run_integral_compare(numpts,True)
