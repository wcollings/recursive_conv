import struct
from sara import Solver

sz_db=8
sz_i=4
def read(fname):
	data=open(fname, 'rb').read()
	offset=0
	order=struct.unpack_from("i",data,0)[0]
	offset+=sz_i
	K_temp=struct.unpack_from("d"*order*2,data,offset)
	K=[]
	for i in range(0,order*2,2):
		K.append(complex(K_temp[i],K_temp[i+1]))
	offset+=sz_db*order*2
	print(K)
	s_temp=struct.unpack_from("d"*order*2,data,offset)
	s=[]
	for i in range(0,order*2,2):
		s.append(complex(s_temp[i],s_temp[i+1]))
	offset+=sz_db*order*2
	print(s)



if __name__=="__main__":
	read("../solv_obj_save.bin")
