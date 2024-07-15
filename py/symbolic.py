from math import factorial, sin,cos,pi,factorial
import numpy as np
from ddiff import ddiff

negsin = lambda x: -sin(x)
negcos = lambda x: -cos(x)
funcs = [ sin, cos, negsin,negcos]
def taylor_coeff(n:int,a:float):
	return funcs[n%4](a)/factorial(n)

def run_ddiff(n:int,a:float):
	step_size = 1/10
	xs = np.linspace(0,n*step_size,n+1)
	ys = np.vectorize(sin)(xs)
	print(xs)
	print(ys)
	derivs = ddiff(xs,ys)
	return derivs

if __name__=="__main__":
	num_terms = 8
	p=pi/2
	tay = 0
	derivs = run_ddiff(num_terms,0)
	print(derivs)
	dd = 0
	print(f"True answer: {sin(p)}")
	for i in range(num_terms):
		tay +=taylor_coeff(i,0)*pow(p,i)
		dd += (derivs[i]/factorial(i))*pow(p,i)
		print(f"{i=}\t{tay=}\t{dd=}")
