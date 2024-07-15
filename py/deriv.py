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

def richardson(f:Callable,x:float, order:int,h:float=0):
	arr = np.zeros((order,1))
	if h==0:
		h=x/2/order
	err = 1e-2
	# Create first set of approximations
	for i in range(order):
		arr[i] = threept(f,x,h/(i+1))
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
	f1 = f(x+2*h)
	f2 = 16*f(x+h)
	f3 = 30*f(x)
	f4 = 16*f(x-h)
	f5 = f(x-2*h)
	return (-f5+f4-f3+f2-f1)/12/h/h

def sevenpt(f,x,h):
	diff = lambda h:f(x+h)-f(x-h)
	e1 = 45*diff(h)
	e2 = 9*diff(2*h)
	e3 = diff(3*h)
	return (e1-e2+e3)/(60*h)

def ninept(f,x,h):
	diff = lambda h:f(x+h)-f(x-h)
	e1 = 864*diff(h)
	e2 = 72*diff(2*h)
	e3 = 32*diff(3*h)
	e4 = 9*diff(4*h)
	return (e1-e2-e3+e4)/1320/h



if __name__=="__main__":
	from poly import Poly
	f=Poly([1,4,5,2,3,5])
	x=3
	o=7
	# print(f"actual deriv: {f.diff()(x)}")
	print(f"Using 3-pt:{threept(L,x,1/10)}")
	print(f"Using 5-pt:{fivept(L,x,1/10)}")
	print(f"Using 7-pt:{sevenpt(L,x,1/10)}")
	print(f"Using 9-pt:{ninept(L,x,1/10)}")
	# print(f"second deriv: {f.diff().diff()(x)}")
	# print(f"using 3-pt:{threept_2(f,x,1/10)}")
	# print(f"Using 5-pt:{fivept_2(f,x,1/10)}")
