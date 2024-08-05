from copy import deepcopy
from cmath import sqrt
from pascal_test import pascal, sumup
def sign(a):
	if a>0:
		return 1
	if a<0:
		return -1
	return 0

def mul_all(l:list[float],c:float):
	return list(map(lambda x:x*c,l))
class Poly:
	"""
	Holds a polynomial. All terms are assumed to be in descending order, 
	i.e. the coefficient of the highest term first, then the second highest, etc.
	"""
	coeff:list[float]
	print_as_tuple:bool
	print_as_you_go:bool
	saved_as_roots:bool
	def __init__(self,c:list,pat:bool=False,pag:bool=False,roots=False):
		self.over_coeff=1.
		self.coeff=[coeff for coeff in c]
		self.print_as_tuple=pat
		self.print_as_you_go=pag
		self.saved_as_roots = roots

	def diff(self):
		out=Poly([0.]*(len(self)-1))
		term=len(self)-1
		for i in range(term):
			out.coeff[i] = self.coeff[i]*(term-i)
		return out
		
	def get_roots(self) -> 'Poly':
		if self.saved_as_roots:
			return self
		if len(self) == 3:
			return bi_find_roots(self)
		elif len(self) == 4:
			return cubic_find_roots_new(self)
		elif len(self)==5:
			return quart_find_roots(self)
		else:
			raise IndexError(f"I don't know how to find the roots of a {len(self)-1} degree polynomial!")

	def depress(self) ->'Poly':
		if self.coeff[0] == 0:
			return deepcopy(self)
		scale = self.coeff[0]
		p=Poly([t/scale for t in self.coeff])
		p.over_coeff=scale
		return p

	def recenter(self,c):
		N = len(self)
		tri = pascal(N)
		dst = [0.] * N
		for i in range(N):
			term=0
			for j in range(i+1):
				loc = sumup(N-j-1)+(i-j)
				bn = self.coeff[j]
				temp=tri[loc]*bn*pow(c,i-j)
				print(tri[loc])
				term += temp
			dst[i] = term
		self.coeff=dst
		return self
	def __mul__(self,rhs:'Poly|float'):
		if isinstance(rhs,float):
			return Poly(mul_all(self.coeff,rhs),self.print_as_tuple)
		num_terms=len(self.coeff)+len(rhs.coeff)-1
		result=[0.0]*num_terms
		for (i,a) in enumerate(self.coeff):
			for (j,b) in enumerate(rhs.coeff):
				result[i+j]+=(a*b)
		# if self.print_as_you_go:
		# 	print(result)
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
		# if self.print_as_you_go:
		# 	print(result)
		return Poly(result,self.print_as_tuple,self.print_as_you_go)
	def __sub__(self,rhs):
		return self + rhs*(-1.0)
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
		# if self.print_as_you_go:
		# 	print(result)
		return Poly(result,self.print_as_tuple,self.print_as_you_go)
	def __repr__(self):
		if self.saved_as_roots:
			res = []
			last_val = None
			for v in self.coeff[::-1]:
				if last_val == v:
					res[-1]+="^2"
				if isinstance(v,complex):
					res.append(f"(x{-v.real:+1.3}{-v.imag:+1.3}j)")
				else:
					res.append(f"(x{-v:+})")
				last_val = v
			return "".join(res)
		if self.print_as_tuple:
			return str(self.coeff)
		res=[]
		temp=""
		pow=len(self.coeff)-1
		temp=f"{self.coeff[0]:-1.3e}"
		if pow > 0:
			temp+="x^"+str(pow)
		res.append(temp)
		for (i,v) in enumerate(self.coeff[1:]):
			temp=""
			pow=len(self.coeff)-i-2
			temp=f"{v:+1.3e}"
			if pow > 0:
				temp+="x^"+str(pow)
			res.append(temp)
		# i=len(self.coeff)-1
		# v=self.coeff[-1]
		# res.append(f"{v}x^{i}")
		return "".join(res)
	def __call__(self,x:float) -> float:
		if self.saved_as_roots:
			res=self.over_coeff
			for i in self.coeff:
				print(f"(x+{-i:-f})={x-i}")
				res *= (x-i)
			return res
		res=self.coeff[0]
		for v in self.coeff[1:]:
			res=v+(x*res)
		return res*self.over_coeff
	def __pow__(self,p):
		if p==1:
			return self
		return self*pow(self,p-1)
	def __len__(self):
		return len(self.coeff)

def bi_find_roots(p:Poly):
	out=Poly(p.coeff)
	scale=out.coeff[0]
	if out.coeff[0] != 1:
		out.coeff=[c/out.coeff[0] for c in out.coeff]
	out.coeff[1]/=2
	z=[]
	b = out.coeff[1]
	c = out.coeff[2]
	if b > 0:
		z.append(-(b+sqrt(b**2-c)))
	else:
		z.append(-b+sqrt(b**2-c))
	z.append(c/z[0])
	p=Poly(z,roots=True)
	p.over_coeff=scale
	return p

def cubic_find_roots_new(p:Poly):
	assert len(p)==4
	out=Poly(p.coeff)
	if out.coeff[0] != 1:
		out.coeff=[c/out.coeff[0] for c in out.coeff]
	shift_factor=0
	if out.coeff[1] != 0:
		b=out.coeff[1]
		shift_factor = b/3
		out.recenter(b/3)
	# cubic is now depressed
	a = out.coeff[2]
	b = out.coeff[3]
	sel = -(b/2)**2 <= (a/3)**3
	print("this eq is ", end='')
	if sel:
		print("easier")
	else:
		print("harder")
	Q = a/3
	R = b/2
	D = Q**3 + R**2
	S = (R + sqrt(D))**(1/3)
	T = (R - sqrt(D))**(1/3)
	print(f'{D=}')
	print(f'{S=}')
	print(f'{T=}')
	z0 = (S+T)
	z1 = complex(-(S+T)/2,sqrt(3)/2*(S-T))
	z2 = complex(-(S+T)/2,-sqrt(3)/2*(S-T))
	# D=complex(-b/2, sqrt(-(a/3)**3-(b/2)**2))
	# z=D**(1./3)
	# print(f"{z=}")
	# z0=2*z.real
	# z1=-z.real + sqrt(3)*z.imag
	# z2=-z.real - sqrt(3)*z.imag
	return Poly([z0,z1,z2],roots=True)
	
def synth_div(p:Poly, z:float):
	res= Poly([])
	last = p.coeff[0]
	res.coeff.append(last)
	for c in p.coeff[1:-1]:
		res.coeff.append(last*z+c)
		last = res.coeff[-1]
	return res

	
def cubic_find_roots(p:Poly):
	assert len(p)==4
	out=Poly(p.coeff)
	if out.coeff[0] != 1:
		out.coeff=[c/out.coeff[0] for c in out.coeff]
	shift_factor=0
	if out.coeff[1] != 0:
		b=out.coeff[1]
		shift_factor = b/3
		out.recenter(b/3)
	bp = -out.coeff[2]
	cp = out.coeff[3]
	start_locs = []
	zeros = []
	print(f"recentered: {out}")
	if bp > 0:
		if abs(cp) == 2*(bp/3)**1.5:
			zeros.append(-sqrt(bp/3)*sign(cp))
			start_locs.append(cp/2/bp + sign(cp)*sqrt(bp))
		else:
			if abs(cp) < 2*(bp/3)**1.5:
				if cp==0:
					start_locs.append(1/2/bp + sqrt(bp))
					start_locs.append(1/2/bp - sqrt(bp))
					start_locs.append(-cp/bp)
				else:
					start_locs.append(cp/2/bp + sign(cp)*sqrt(bp))
					start_locs.append(cp/2/bp - sign(cp)*sqrt(bp))
					start_locs.append(-cp/bp)
			else:
				if cp**2 > abs(bp**3):
					start_locs.append(cp**(1/3))
				else:
					start_locs.append(cp/2/bp + sign(cp)*sqrt(bp))
	else:
		if cp**2 > abs(bp**3):
			start_locs.append(cp**(1/3))
		else:
			start_locs.append(-cp/2/bp)
	if zeros:
		print(f"Found a zero early: {zeros[0]}")
	for pt in start_locs:
		print(f"looking around {pt}")
		zeros.append(newton(p,pt) -shift_factor)
	return Poly(zeros,roots=True)

def quart_find_roots(p:Poly):
	p.print_as_tuple = True
	assert len(p)==5
	shift_factor=0
	out = Poly(p.coeff,pat=True)
	if out.coeff[0] != 1:
		out.coeff=[c/out.coeff[0] for c in out.coeff]
	if out.coeff[1] != 0:
		shift_factor = out.coeff[1]/4
		out.recenter(out.coeff[1]/4)
	c=out.coeff[2]
	d=out.coeff[3]
	e=out.coeff[4]
	print("shifted:")
	print(out)
	tri = Poly([1,2*c,c**2-4*e,-d**2], pat=True)
	print("tri:")
	print(tri)
	print("Solved:")
	print(cubic_find_roots(tri))
	z0 = sqrt(abs(cubic_find_roots(tri).coeff[0]))
	z1 = -z0
	z2 = (c+z0**2 + d/z0)/2
	z3 = (c+z0**2 - d/z0)/2
	roots = [z0,z1,z2,z3]
	return Poly([root+shift_factor for root in roots],roots=True)
		

def newton(f:Poly,xn:float,err:float=1e-8):
	fp = f.diff()
	while abs(f(xn)) > err:
		xn = xn - (f(xn)/fp(xn))
	return xn

if __name__=="__main__":
	d=Poly([5,11,1])
	d1=d.get_roots()
	d1.coeff=[s.real for s in d1.coeff]
	print(d(4))
	print(d1(4))
	print(synth_div(d,d1.coeff[0]))
