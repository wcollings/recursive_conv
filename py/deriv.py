from typing import Callable
from ddiff import L
import numpy as np
def richardson_lut(arr:np.ndarray):
	order = arr.shape[0]
	err = 1e-2
	final_idx=order-1
	for i in range(1,order):
		for j in range(order-1,i-1,-1):
			denom = max(pow(2,i-1)-1,1)
			arr[j] = arr[j] + (arr[j]-arr[j-1])/denom
		if abs(arr[i]-arr[i-1]) < err:
			final_idx = i
			print("Finished early!")
			break
	return arr[final_idx]

def richardson(f:Callable,x:float, order:int,h:float=1/10):
	arr = np.zeros((order,1))
	if h==0:
		h=x/2/order
	err = 1e-2
	# Create first set of approximations
	for i in range(order):
		arr[i] = ninept_7(f,x,h/(i+1))
	# Refine
	final_idx = order-1
	for i in range(1,order):
		for j in range(order-1,i-1,-1):
			denom = max(pow(2,i-1)-1,1)
			arr[j] = arr[j] + (arr[j]-arr[j-1])/denom
		if abs(arr[i]-arr[i-1]) < err:
			final_idx = i
			print("Finished early!")
			break
	return arr[final_idx]

def threept(f,x,h):
	return (f(x+h)-f(x-h))/(2*h)
def threept_2(f,x,h):
	return (f(x+2*h)-2*f(x)+f(x-2*h))/h/h/4

def fivept(f,x,h):
	diff = lambda h:f(x+h)-f(x-h)
	e1 = 8*diff(h)
	e2 = diff(2*h)
	return (e1-e2)/12/h

def fivept_2(f,x,h):
	diff = lambda h:f(x+h)+f(x-h)
	e0 = 30*f(x)
	e1 = 16*diff(h)
	e2 = diff(2*h)
	return (e1-e2-e0)/12/h/h

def sevenpt(f,x,h):
	diff = lambda h:f(x+h)-f(x-h)
	e1 = 45*diff(h)
	e2 = 9*diff(2*h)
	e3 = diff(3*h)
	return (e1-e2+e3)/(60*h)

def sevenpt_2(f,x,h):
	s = lambda h:f(x+h)+f(x-h)
	e0=120*f(x)
	e1=25*s(h)
	e2=44*s(2*h)
	e3=9*s(h*3)
	return (e1+e2-e3-e0)/120/h/h

def sevenpt_3(f,x,h):
	diff = lambda h:f(x+h)-f(x-h)
	e1=13*diff(h)
	e2=8*diff(2*h)
	e3=diff(3*h)
	return (-e1+e2-e3)/8/h/h/h


def ninept(f,x,h):
	diff = lambda h:f(x+h)-f(x-h)
	e1=1050*diff(h)
	e2=300*diff(2*h)
	e3=75*diff(3*h)
	e4=25*diff(4*h)/2
	e5=diff(5*h)
	return (e1-e2+e3-e4+e5)/1260/h

def ninept_2(f,x,h):
	s = lambda h:f(x+h)+f(x-h)
	e0=45360*f(x)
	e1=5502*s(h)
	e2=61656*s(2*h)
	e3=27621*s(3*h)
	e4=6488*s(4*h)
	e5=665*s(5*h)
	return (e1+e2-e3+e4-e5-e0)/90720/h/h

def ninept_3(f,x,h):
	diff = lambda h:f(x+h)-f(x-h)
	coeffs = [70089,-52428,14607,-2522,-205]
	tot = sum(c * diff(h*i) for i,c in enumerate(coeffs))
	return tot/-30240/h/h/h
def ninept_4(f,x,h):
	s = lambda h:f(x+h)+f(x-h)
	pass
def ninept_5(f,x,h):
	s = lambda h:f(x+h)+f(x-h)
	pass
def ninept_6(f,x,h):
	s = lambda h:f(x+h)+f(x-h)
	pass
def ninept_7(f,x,h):
	diff = lambda h:f(x+h)-f(x-h)
	e1=378*diff(h)
	e2=408*diff(2*h)
	e3=207*diff(3*h)
	e4=52*diff(4*h)
	e5=5*diff(5*h)
	return (e1-e2+e3-e4+e5)/-24/pow(h,7)

def eval_stencil(f,x,h,p,coeffs):
	tot = 0
	denom = coeffs.pop()
	if p%2:
		sd = lambda h:f(x+h)-f(x-h)
	else:
		sd=lambda h:f(x+h)+f(x-h)
		tot=coeffs[0]*f(x)
		coeffs.pop(0)
	tot += sum(c*sd(h*(i+1)) for i,c in enumerate(coeffs))
	# breakpoint()
	return tot/(denom*pow(h,p))

stencils = {
	3:[
		[1,2],
		],
	5:[
		[8,-1,12],
		[-30,16,-1,12]
		],
	7:[
		[35,-56,28,-8,1,-70],
		[-120,25,44,-9,120],
		[-13,8,-1,8],
		[0,899,-1612,837,-124,930],
		[5,-4,1,2],
		[0,-217,434,-279,62,217],
		[-14,14,-6,1,2]
		],
	9:[]
}

if __name__=="__main__":
	from poly import Poly
	f=Poly([1,4,5,2,3,5,4,2])
	ordinals = ["first","second","third","fourth","fifth","sixth","seventh"]
	x=3
	o=7
	h=1/10
	deriv=Poly([1,4,5,2,3,5,4,2])
	for i in range(7):
		deriv=deriv.diff()
		print(f"{ordinals[i]} derivative")
		print(f"true answer: {deriv(x)}")
		for j in [3,5,7,9]:
			if len(stencils[j]) > i:
				approx=eval_stencil(f,x,h,i+1,stencils[j][i])
				err = abs((deriv(x)-approx))/approx
				print(f"using {j}-pt:  {approx:.4f} ({err=:1.2e})")
		print("-------------")
	print("Derivatives of L:")
	for i in range(7):
		j=7
		approx=eval_stencil(L,x,h,i+1,stencils[j][i])
		print(f"{ordinals[i]}:\t{approx:.4e}")
