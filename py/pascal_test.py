from math import floor
def sumup(row):
	return int(row*(row+1)/2)

def print_arr(arr,arrows):
	for v in arr:
		print(f"{v: <3}",end="")
	print()
	idxs=['   ']*len(arr)
	for arrow in arrows:
		idxs[arrow]='^  '
	print("".join(idxs))
	
def pascal(row):
	ne = sumup(row+1)
	# print(ne)
	tri = [0]*ne
	col=0
	for i in range(ne):
		start_col = sumup(col)
		end_col = sumup(col+1)-1
		# print(f"{col=}",end=" ")
		if (i in [start_col,end_col]):
			tri[i]=1
			if i == end_col:
				# print(f"row={end_col}")
				col+=1
			# else:
				# print(f"row=0 BOUNDARY")
			# print_arr(tri,[i])
		else:
			offset=i-sumup(col)
			# print(f"row={offset}")
			idx=sumup(col-1)+offset-1
			tri[i]=tri[idx]+tri[idx+1]
			# print(f"{tri[idx]} + {tri[idx+1]} ({idx=})")
			# print_arr(tri,[i,idx,idx+1])
	for i in range(1,row+1):
		for j in range(1,floor((i-1)/2)+2):
			tri[sumup(i)+2*j-1]*= -1
	# print(tri)
	return tri
def recenter(l,c):
	N = len(l)
	tri = pascal(N)
	dst = [0] * N
	for i in range(N):
		term=0
		for j in range(i+1):
			loc = sumup(N-j-1)+(i-j)
			bn = l[j]
			temp=tri[loc]*bn*pow(c,i-j)
			term += temp
		dst[i] = term

if __name__=="__main__":
	recenter([1,2,5],3)
