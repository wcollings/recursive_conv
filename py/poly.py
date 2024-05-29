from copy import deepcopy
def mul_all(l:list[float],c:float):
	return list(map(lambda x:x*c,l))
class Poly:
	coeff:list[float]
	print_as_tuple:bool
	print_as_you_go:bool
	def __init__(self,c:list,pat:bool=False,pag:bool=False):
		self.coeff=c
		self.print_as_tuple=pat
		self.print_as_you_go=pag
	def __mul__(self,rhs:'Poly|float'):
		if isinstance(rhs,float):
			return Poly(mul_all(self.coeff,rhs),self.print_as_tuple)
		num_terms=len(self.coeff)+len(rhs.coeff)-1
		result=[0.0]*num_terms
		for (i,a) in enumerate(self.coeff):
			for (j,b) in enumerate(rhs.coeff):
				result[i+j]+=(a*b)
		if self.print_as_you_go:
			print(result)
		return Poly(result,self.print_as_tuple,self.print_as_you_go)
	def __add__(self,rhs:'Poly'):
		better=max(self.coeff,rhs.coeff,key=len)
		left=(better==self.coeff)
		num_terms=len(better)
		if left:
			result=deepcopy(self.coeff)
			to_add=rhs.coeff
		else:
			result=deepcopy(rhs.coeff)
			to_add=self.coeff
		for i,v in enumerate(to_add):
			result[i]+=v
		if self.print_as_you_go:
			print(result)
		return Poly(result,self.print_as_tuple,self.print_as_you_go)
	def __radd__(self,rhs:'Poly'):
		better=max(self.coeff,rhs.coeff,key=len)
		left=(better==self.coeff)
		num_terms=len(better)
		if left:
			result=deepcopy(self.coeff)
			to_add=rhs.coeff
		else:
			result=deepcopy(rhs.coeff)
			to_add=self.coeff
		for i,v in enumerate(to_add):
			result[i]+=v
		if self.print_as_you_go:
			print(result)
		return Poly(result,self.print_as_tuple,self.print_as_you_go)
	def __repr__(self):
		if self.print_as_tuple:
			return str(self.coeff)
		res=[]
		for (i,v) in enumerate(self.coeff):
			temp=""
			pow=len(self.coeff)-i-1
			temp=f"{v:+}"
			if pow > 0:
				temp+="x^"+str(pow)
			res.append(temp)
		# i=len(self.coeff)-1
		# v=self.coeff[-1]
		# res.append(f"{v}x^{i}")
		return "".join(res)
	def __call__(self,x:float) -> float:
		res=self.coeff[0]
		for v in self.coeff[1:]:
			# print(f"{v}x^{i}={v*(x**i)}")
			res=v+(x*res)
		return res
	def __pow__(self,p):
		if p==1:
			return self
		return self*pow(self,p-1)
	def __len__(self):
		return len(self.coeff)-1

if __name__=="__main__":
	p1=Poly([1,-1])
	p2=Poly([1,-2])
	p3=Poly([1,-3])
	p4=Poly([1,-4])
	p12=p1*p2
	p122=p12*p2
	p1223=p122*p3
	p12233=p1223*p3
	p122334=p12233*p4
	p0=Poly([2])
	print(p1.coeff)
	print(p1)
	p0.print_as_you_go=False
	p12=p1*p2
	p122=p1*p2*p2*-2.0
	p1223=p1*p2*p2*p3*1.5
	p12233=p1*p2*p2*p3*p3*-1.5
	p122334=p1*p2*p2*p3*p3*p4*0.75
	print(p0+p12+p122+p1223+p12233+p122334)
