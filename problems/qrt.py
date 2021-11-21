with open('rosalind_qrt.txt') as f:
	lines = f.readlines()

taxa = lines[0].strip().split(' ')
quartets = set()
for line in lines[1:]:
	C = []
	D = []
	for i in range(len(line)):
		if line[i] == '1':
			C.append(taxa[i])
		elif line[i] == '0':
			D.append(taxa[i])
	if len(C) >= 2 and len(D) >= 2:
		for i in range(len(C)-1):
			for j in range(i+1, len(C)):
				for k in range(len(D)-1):
					for l in range(k+1, len(D)):
						if C[i] < C[j]:
							A = (C[i], C[j])
						else:
							A = (C[j], C[i])
						if D[k] < D[l]:
							B = (D[k], D[l])
						else:
							B = (D[l], D[k])
						quartets.add((A, B) if A[0] < B[0] else (B, A))

for quartet in quartets:
	A, B = quartet
	print('{{{}, {}}} {{{}, {}}}'.format(A[0], A[1], B[0], B[1]))