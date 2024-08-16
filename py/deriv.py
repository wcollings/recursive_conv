from typing import Callable
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
	e2=62656*s(2*h)
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
	denom = coeffs[-1]
	coeff = coeffs[:-1]
	# denom = coeffs.pop()
	if p%2:
		sd = lambda h:f(x+h)-f(x-h)
	else:
		sd=lambda h:f(x+h)+f(x-h)
		tot=coeffs[0]*f(x)
		coeff=coeffs[1:-1]
	tot += sum(c*sd(h*(i+1)) for i,c in enumerate(coeff))
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
		# [35,-56,28,-8,1,-70],
		[45,-9,1,60],
		[-120,25,44,-9,120],
		[-13,8,-1,8],
		[0,-13,16,-3,240],
		# [0,899,-1612,837,-124,930],
		[5,-4,1,2],
		[0,-217,434,-279,62,217],
		[-14,14,-6,1,2]
		],
	9:[
		[864,-72,-32,9,1320],
		[-45360,5502,62656,-27621,6488,-665,90720],
		[70089,-52428,14607,-2522,205,-30240],
		[0,7322,-10336,6081,-1552,165,3360],

		]
}

# L_derivs = list(map(float,open("L_deriv.csv").readline().split(',')))

def take_derivs(l:Callable,x0:float,num:int,h:float):
	res = [l(x0)]
	for d in range(num):
		res.append(eval_stencil(l,x0,h,d+1,stencils[7][d]))
	return res

if __name__=="__main__":
	from poly import Poly
	c = Poly(L_derivs[::-1])
	ordinals = ["first","second","third","fourth","fifth","sixth","seventh"]
	# print(sevenpt(L,2e6,1/10))
	x0=1e8
	x1=2e5
	o=7
	h=10
	l = lambda s:L(s)-L(x0)
	approx=eval_stencil(l,x1,2,3,stencils[7][0])
	print(f"zeroth:\t{L(x1)-L(x0):.4e}")
	for i in range(7):
		approx=eval_stencil(l,x1,h,i+1,stencils[7][i])
		err = abs((L_derivs[i]-approx)/L_derivs[i])
		print(f"{ordinals[i]}:\t{approx}\t({err:1.4e})")
