from math import exp
K=[0.9478,0.0522]
s=[-2.105,-0.095]
T=40
delta_n=1
y=[]
x=[]
R2=2

def Phi(i):
	return exp(s[i]*delta_n)
def zeta(i):
	return -(s[i]*delta_n)
def q1(i:int):
	q0=(delta_n/zeta(i))*(1-Phi(i))
	return q0,
def q2(i:int):
	zi=zeta(i)
	q0=(delta_n/zi**2)*(-1+zi+Phi(i))
	q1=(delta_n/zi**2)*(1-(1+zi)*Phi(i))
	return q0,q1
def q3(i:int):
	zi=zeta(i)
	q0=(delta_n/(2*zi**3))*(2-(3*zi)+(2*zi**2) - (2-zi)*Phi(i))
	q1=(delta_n/zi**3)*(-2*(1-zi)+(2-zi**2)*Phi(i))
	q2=(delta_n/(2*zi**3))*(2-zi-(2+zi)*Phi(i))
	return q0,q1,q2
qs=[q1,q2,q3]
calc_q=qs[R2-1]
