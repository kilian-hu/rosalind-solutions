with open("rosalind_revc.txt") as f:
	s = f.readline().strip()

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
for x in reversed(range(len(s))):
	print(complements[s[x]], end='')
print()