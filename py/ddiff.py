from typing import Iterable
from math import atan, factorial

def ddiff(x:Iterable, y:Iterable,ne:int=-1):
	out=[i for i in y]
	ne=len(out)
	x=tuple(x)
	for i in range(1,ne+1):
		for j in range(ne-1,i-1,-1):
			# print(f"(y[{j}]-y[{j-1}])/(x[{j}]-x[{j-i}])",end=" = ")
			# print(f"({out[j]:1.4f}-{out[j-1]:1.4f})/({x[j]:1.1f}-{x[j-i]:1.1f})",end =" = ")
			num = out[j]-out[j-1]
			denom = x[j]-x[j-i]
			# print(f"{num:1.4f}/{denom:1.1f}={num/denom:1.4f}")
			out[j]=num/denom
	out[-1]=(out[-1]-out[-2])/(x[-1]-x[0])
	return out

def scale_to_taylor(terms:list[float]):
	return list(map(lambda i:terms[i]/factorial(i),range(len(terms))))

def L(f):
	a=1e-9
	b=2.8e-9
	c=800e-9
	f0=2e4
	res=(0.6366*a)*atan(-c*(f-f0))+b;
	return res
if __name__=="__main__":
	start=int(2e6)
	step=20
	ne=10
	xs=list(range(start,start+(ne*step),step))
	ys=tuple(map(L,xs))
	dd=ddiff(xs,ys)
	print("x\t\ty\t\tdd")
	for a,b,c in zip(xs,ys,dd):
		print(f"{a:e}\t{b:e}\t{c:e}")
	# print(ddiff(xs,ys))
