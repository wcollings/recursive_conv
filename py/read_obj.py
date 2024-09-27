from struct import unpack_from
from dataclasses import dataclass
from typing import Iterable
sz_i=4
sz_db=8
@dataclass
class Solver:
	order:int
	K0:complex
	Ki:Iterable[complex]
	si:Iterable[complex]
	curr_t:float
	curr_x:float
	tt:Iterable[float]
	xx:Iterable[float]
	yy:Iterable[complex]

def read_f(fname):
	contents=open(fname,"rb").read()
	offset=0
	order=unpack_from("i",contents,offset)[0]
	offset+=sz_i
	K0=unpack_from("dd",contents,offset)[0]
	offset+=sz_db*2
	Ki=[]
	for i in range(order):
		Ki.append(complex(*unpack_from("dd",contents,offset)))
		offset+=2*sz_db
	si=[]
	for i in range(order):
		si.append(complex(*unpack_from("dd",contents,offset)))
		offset+=2*sz_db
	t=unpack_from("d",contents,offset)[0]
	offset+=sz_db
	x=unpack_from("d",contents,offset)[0]
	offset+=sz_db
	tt=unpack_from("d"*order,contents,offset)
	offset+=sz_db*order
	xx=unpack_from("d"*order,contents,offset)
	offset+=sz_db*order
	yy=[]
	for i in range(order):
		yy.append(complex(*unpack_from("dd",contents,offset)))
		offset+=2*sz_db
	solv=Solver(order,K0,Ki,si,t,x,tt,xx,yy)
	# Make it print whatever you need

if __name__=="__main__":
	read_f("../solv_output.bin")
