with open("rosalind_tree.txt") as file:
	# Node indexing starts at 1.
	n = int(file.readline().strip())
	adj_list = [[] for _ in range(n+1)]
	for u, v in [map(int, line.strip().split(' ')) for line in file.readlines()]:
		adj_list[u].append(v)
		adj_list[v].append(u)

marked = [False for _ in range(n+1)]
num_components = 0
for i in range(1, n+1):
	if marked[i]:
		continue
	marked[i] = True
	num_components += 1
	open_front = [i]
	while len(open_front) > 0:
		u = open_front.pop()
		unmarked_neighbours = [v for v in adj_list[u] if not marked[v]]
		open_front += unmarked_neighbours
		for v in unmarked_neighbours:
			marked[v] = True

# All components can be connected with num_components - 1 edges.
print(num_components - 1)