import numpy as np
from plotting import figure_wrapper
from matplotlib import pyplot as plt
import pandas as pd
from poly import Poly
from pade import Pade, separate
from typing import Callable
from scipy.optimize import curve_fit as cf

df=pd.read_csv("results/zbus2.csv")
# df['l']*=1e-9
# df['r']*=1e-3
# df['r']+=min(df['r'])
df['w']=df.f*2*np.pi
freq=df.w
df['k0']=1/df.l
df['s0']=df.r/df.l
df['y']=df.k0/(freq+df.s0)
# lower_idx=0
lower_idx=178

def create_pade(n:int,d:int) -> Callable:
	def inner(*vals):
		assert len(vals)==n+d,f"{n=},{d=},{vals=}"
		num=Poly(vals[:n])
		denom=Poly(vals[n:])
		return Pade(num,denom)
	return inner
def call_pade_param(n,d):
	p=create_pade(n,d)
	return lambda f,*a:p(*a)(f)
def run_optim(n,d):
	base = call_pade_param(n,d)
	coeff=cf(base,freq[lower_idx:],df.y[lower_idx:],p0=[1]*(n+d))
	model=create_pade(n,d)(*coeff[0])
	model.ind_var='s'
	# print(f"h{n}{d}:")
	# print(separate(model).get_coeff_str())
	return model
	return df.f.map(model)

h12 = run_optim(1,2)
h23 = run_optim(2,3)
freqs=[10**i for i in range(3,9)]

if __name__=="__main__":
	with open("ybus.txt",'w') as fp:
		fp.write(separate(h23).get_coeff_str())
	Y=df.k0/(df.f+df.s0)
	plt.semilogx(freq,df.y,label="Model")
	plt.semilogx(freq,h23(freq),label="Experimenal")
	plt.xlabel("Frequency (rad/sec)")
	plt.ylabel("Admittance (S)")
	plt.legend()
	plt.show()
	# with figure_wrapper(outf="Ybus.png",interactive=True) as fw:
	# 	fw.slogx(freq,df.y,name="Experiment")
	# 	fw.slogx(freq,h23(freq),name="Model")
	# 	fw.set_labels("frequency (Hz)", "Admittance (Siemans)")
	# 	fw.set_xlim(1e3,1.2e8)
	# 	fw.set_ylim(-5,140)
	# 	fw.set_fontsize(15)
	# 	fw.make_legend=True
	# 	# fw.ylabel("Admittance (Mho)")
		# fw.grid(True)
