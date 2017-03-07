b = [0,6,1,7,4]

def  maxDifference(b):
	b = sorted(set(b), key=b.index)
	a = list(reversed(b))
	highest = 0
	for num in range(len(b)):
		if num != (len(b)-1):
			for comparenum in range(num+1,len(b)):
				if a[num] - a[comparenum] > highest:
					highest = a[num] - a[comparenum]
	if highest == 0:
		return -1
	else:
		return highest
		


print maxDifference([2,3,10,2,4,8,1])