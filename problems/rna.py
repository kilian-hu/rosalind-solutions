with open("rosalind_rna.txt") as f:
	t = f.readline().strip()

print(t.replace('T', 'U'))