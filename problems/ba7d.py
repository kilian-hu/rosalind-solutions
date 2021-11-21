import sys


class Node():
	def __init__(self, nid, age, size, p=None, l=None, r=None):
		self.nid = nid  # Node ID.
		self.age = age  # Age of the node.
		self.size = size  # Cluster size.
		self.p = p  # Parent.
		self.l = l  # Left child.
		self.r = r  # Right child.

	def __repr__(self):
		return "nid={}\tage={}\tsize={}\tp={}\tl={}\tr={}".format(self.nid, self.age, self.size, self.p, self.l, self.r)


with open("rosalind_ba7d.txt") as f:
	n = int(f.readline().strip())
	D = dict()
	for i in range(n):
		D[i] = dict(zip(range(n), [float(x) for x in f.readline().strip().split('\t') if x != '']))

G = [Node(nid=i, age=0.0, size=1) for i in range(n)]

while len(D) > 1:
	# Find smallest distance.
	best_d, best_i, best_j = sys.float_info.max, None, None
	for i, i_row in D.items():
		for j, d in i_row.items():
			if d < best_d and i != j:
				best_d, best_i, best_j = d, i, j
	if best_i > best_j:
		best_i, best_j = best_j, best_i
	new_nid = len(G)
	new_size = G[best_i].size + G[best_j].size
	G.append(Node(nid=new_nid, age=best_d/2, size=new_size, l=best_i, r=best_j))
	G[best_i].p = new_nid
	G[best_j].p = new_nid
	new_row = {new_nid: 0.0}
	i_row = D[best_i]
	j_row = D[best_j]
	del D[best_i]
	del D[best_j]
	for k, k_row in D.items():
		new_d = (G[best_i].size*i_row[k] + G[best_j].size*j_row[k]) / G[new_nid].size
		new_row[k] = new_d
		D[k][new_nid] = new_d
		del k_row[best_i]
		del k_row[best_j]
	D[new_nid] = new_row

for node in G:
	if node.l != None:
		print("{}->{}:{:.3f}".format(node.nid, node.l, node.age-G[node.l].age))
		print("{}->{}:{:.3f}".format(node.nid, node.r, node.age-G[node.r].age))
	if node.p != None:
		print("{}->{}:{:.3f}".format(node.nid, node.p, G[node.p].age-node.age))