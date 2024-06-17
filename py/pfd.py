def pow(a:float,b:float) -> float:
	return a**b
def eval_p(t,z) -> float:
	terms=[6,11,6,1]
	res=0
	for i,term in enumerate(terms):
		res+=term*pow(t,i)
	return res/(t-z)

def eval_sd(t,z):
	terms=[1,6,11,6]
	div_terms=[1]
	for term in terms[1:]:
		div_terms.append(div_terms[-1]*z+term)
	res=0
	div_terms=div_terms[::-1][1:]
	print(div_terms)
	for i,term in enumerate(div_terms):
		res+=term*pow(t,i)
	return res

print(eval_p(14,-3))
print(eval_sd(14,-3))
