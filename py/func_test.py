from math import tanh
import pyqtgraph as pg
from numpy import linspace
def fac(x:int):
	if x <0:
		raise ValueError("Negative values are not allowed!")
	if x==0:
		return 1
	return x*fac(x-1)
def sinh(x:float,nt:int):
	tot = 0
	for n in range(nt):
		p=2*n+1
		tot += (x**p)/fac(p)
	return tot

def cosh(x:float,nt:int):
	tot = 0
	for n in range(nt):
		p=2*n
		tot += (x**p)/fac(p)
	return tot

if __name__=="__main__":
	true_vals = []
	xs = linspace(-10,10,50)
	nt = 4
	approx = []
	err = []
	for x in xs:
		a = tanh(-x)
		b = sinh(-x,nt-1)/cosh(-x,nt)
		true_vals.append(a)
		approx.append(b)
		err.append(abs(a-b)/a)
	app = pg.mkQApp()
	win = pg.GraphicsLayoutWidget(show=True)
	p1 = win.addPlot()
	p1.plot(xs,true_vals, pen = (255,0,0))
	p1.plot(xs,approx, pen = (0,255,0))
	p1.showGrid(True,True)
	win.nextRow()
	p2 = win.addPlot()
	p2.plot(xs,err)
	p2.showGrid(True,True)
	app.exec()
