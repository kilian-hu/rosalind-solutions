with open('rosalind_cntq.txt') as f:
	lines = f.readlines()

n = int(lines[0].strip())
memo = [1] + [None for _ in range(n)]
for i in range(1, n+1):
	memo[i] = i * memo[i-1]
print((memo[n]/(memo[n-4]*memo[4]))%1000000)