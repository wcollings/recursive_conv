from matplotlib import pyplot as plt
from poly import Poly
import numpy as np
import pandas as pd

class quad:
	def __init__(self) -> None:
		self.val=0
		self.last=0
		self.last_time=0
	def __call__(self,t,v):
		_last = self.val
		self.val+=(t-self.last_time)*(v+self.last)/2
		self.last=_last
		self.last_time=t
		return self.val
	
def conv(x,y):
	res = []
	for i in range(len(y)):
		s=0
		for j in range(min(len(x),i)):
			top=len(x)-1
			s+=y[i-j]*x[top-j]
		res.append(s)
	return res

if __name__=="__main__":
	xs=np.arange(0,40,1e-2)
	p=Poly([1,2,2,1])
	dp=p.diff()
	int_v = []
	approx=quad()
	res=[]
	for x in xs:
		y=p(x)
		yp=dp(x)
		res.append(approx(x,yp))
		int_v.append(y)
	plt.plot(xs,int_v)
	plt.plot(xs,res)
	plt.show()
