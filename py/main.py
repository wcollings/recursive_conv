from deriv import der
from pade import solve_system,eval_pade
from poly import Poly
from math import atan
import numpy as np
from typing import Callable
import pyqtgraph as pg

def L(f,f0=0):
	a=1e-9
	b=2.8e-9
	c=800e-9
	res=(0.6366*a)*atan(-c*(f-f0))+b;
	return res

def eval_pade_p(num:Poly,denom:Poly,x:float):
	return num(x)/denom(x)
def test_fit(exp:Callable,approx:Callable, xs:np.ndarray):
	exp=np.vectorize(exp)
	ys_e = exp(xs)
	ys_a = approx(xs)
	err = abs(ys_e-ys_a)/ys_e
	approx=np.vectorize(approx)
	app = pg.mkQApp()
	win = pg.GraphicsLayoutWidget(show=True)
	p1 = win.addPlot()
	p1.plot(xs,ys_e, pen = (255,0,0))
	p1.plot(xs,ys_a, pen = (0,255,0))
	p1.setLogMode(True,True)
	p1.showGrid(True,True)
	win.nextRow()
	p2 = win.addPlot()
	p2.plot(xs,err)
	p2.setLogMode(True,True)
	p2.showGrid(True,True)
	app.exec()

def plot_L():
	H1 = np.vectorize(lambda s: Poly([5,1])(s)/Poly([5,11,1])(s))
	app = pg.mkQApp()
	win = pg.GraphicsLayoutWidget(show=True)
	p1 = win.addPlot()
	xs = np.logspace(1,5)
	ys = H1(xs)
	p1.plot(xs,ys, pen = (255,0,0))
	p1.setLogMode(True,True)
	p1.showGrid(True,True)
	app.exec()

if __name__=="__main__":
	start=int(2e6)
	step=100
	ne=10
	xs=list(range(start,start+(ne*step),step))
	for i in range(7):
		j=7
		approx=eval_stencil(L,x,h,i+1,stencils[j][i])
	L1 = lambda x:L(x,start)
	ys=tuple(map(L1,xs))
	dd=scale_to_taylor(ddiff(xs,ys))
	# for a,b,c in zip(xs,ys,dd):
	# 	print(f"{a:e}\t{b:e}\t{c:e}")
	aa,bb=solve_system(dd,3,4)
	num = Poly(aa[::-1]) #pyright:ignore
	denom = Poly(bb[::-1]) #pyright:ignore
	# taylor = Poly(dd[::-1]) #.recenter(start)
	test_xs = np.logspace(3,8,1000)
	approx = lambda x:eval_pade_p(num,denom,x)
	test_fit(L1,approx,test_xs)
