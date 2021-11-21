from collections import defaultdict


with open("rosalind_ba7a.txt") as f:
	n = int(f.readline().strip())
	adj = defaultdict(list)
	lines = f.readlines()
	for line in lines:  # u->v:w
		if line == '\n':
			continue
		u, tmp = line.strip().split("->")  # [u, v:w]
		v, w = tmp.split(":")
		adj[int(u)].append((int(v), int(w)))

for i in range(n):
	i_res = [None for j in range(n)]
	i_res[i] = 0
	marked = dict(zip(adj.keys(), [False for _ in range(len(adj))]))
	marked[i] = True
	frontier = [vw for vw in adj[i]]
	while len(frontier) > 0:
		v, w = frontier.pop(0)
		marked[v] = True
		if v < n:
			i_res[v] = w
		frontier += [(u, w+x) for u, x in adj[v] if not marked[u]]
	for r in i_res[:-1]:
		print(r, end=' ')
	print(i_res[-1])