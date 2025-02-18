import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from poly import Poly
from pade import Pade, separate
from typing import Callable
from scipy.optimize import curve_fit as cf
from plotting import figure_wrapper

df=pd.read_csv("results/zbus2.csv")
ladder=pd.read_csv("results/Y_ladder.csv")
ladder['w']=ladder.f*2*np.pi
# df['l']*=1e-9
# df['r']*=1e-3
# df['r']+=min(df['r'])
df['w']=df.f*2*np.pi
freq=df.w
df['k0']=1/df.l
df['s0']=df.r/df.l
df['y']=df.k0/(freq+df.s0)
lower_idx=0
# lower_idx=178

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
	with figure_wrapper(interactive=True,outf="Y_vs_ladder.png") as fw:
		fw.slogx(freq,df.y,name="Experiment",color="black")
		fw.slogx(ladder.w,ladder.y,name="Ladder Network",color="#0072BD")
		fw.slogx(freq,h23(freq),name="Model",linestyle="--",color="#cc0000")
		fw.xlim=1e4,1e9
		fw.ylim=-7,166
		fw.set_fontsize(20)
		fw.set_labels("Frequency (rad/sec)","|Admittance| (S)")

	# plt.semilogx(freq,df.y,label="Model")
	# plt.semilogx(ladder.w,ladder.y,label="Ladder Network")
	# plt.semilogx(freq,h23(freq),label="Experimenal")
	# plt.xlabel("Frequency (rad/sec)")
	# plt.ylabel("Admittance (S)")
	# plt.grid()
	# plt.legend()
	# plt.show()
