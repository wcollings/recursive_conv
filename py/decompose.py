from math import sqrt,cbrt
from poly import Poly
def roots_3(p:Poly):
	assert len(p)==3,"You didn't give me a cubic polynomial!"
	d=p.coeff[0]
	c=p.coeff[1]
	b=p.coeff[2]
	a=p.coeff[3]
	delta_0=b**2-(3*a*c)
	delta_1=2*b**3 - 9*a*b*c - 27*a**2*d
	e=(-1+complex(0,sqrt(3)))/2
	C=cbrt((delta_1+sqrt(delta_1**2-(4*delta_0**3)))/2)
