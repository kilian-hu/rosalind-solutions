with open("rosalind_dna.txt") as f:
	s = f.readline().strip()

counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
for x in s:
	counts[x] += 1

print(counts['A'], counts['C'], counts['G'], counts['T'])