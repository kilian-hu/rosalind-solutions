VALID_CHARS = {'A', 'T', 'C', 'G'}

class DNA():
	def __init__(self, name='', info='', content=None):
		self.name = name
		self.info = info
		self.content = content

	def __repr__(self):
		return "name: {}\ninfo: {}\ncontent: {}".format(self.name, self.info, self.content)


def add_new_DNA(dna_list, line):
	assert line[0] == '>'
	first_space_idx = line.find(' ')
	if first_space_idx != -1:
		dna_name = line[1:first_space_idx]
		dna_info = line[first_space_idx:].strip()
	else:
		dna_name = line[1:]
		dna_info = ''
	dna_list.append(DNA(name=dna_name, info=dna_info, content=[]))


def add_line_to_DNA(cur_DNA, line):
	for x in line:
		if x in VALID_CHARS:
			cur_DNA.content.append(x)
		elif x == ' ':
			continue
		else:
			raise Exception()


def parse_FASTA(file):
	"""
	Basic state machine for parsing.
	0 = Expecting '>' or empty line.
	1 = Expecting valid string for current DNA or empty line.
	2 = Got at least one string for current DNA. Ready for new DNA
	or continue current DNA.
	"""
	state = 0
	dna_list = []
	for line in file:
		line = line.strip()
		if state == 0:
			if line[0] == '>':
				add_new_DNA(dna_list, line)
				state = 1
			elif line == '':
				continue
			else:
				raise Exception()
		elif state == 1:
			add_line_to_DNA(dna_list[-1], line)
			state = 2
		elif state == 2:
			if line[0] == '>':
				add_new_DNA(dna_list, line)
				state = 1
			else:
				add_line_to_DNA(dna_list[-1], line)
		else:
			raise Exception()
	file.seek(0)
	return dna_list


# GO HERE
with open('rosalind_mult.txt') as file:
	dna_list = parse_FASTA(file)

a, b, c, d = [x.content for x in dna_list]

def get_best_score(D, i, j, k, l, ax, bx, cx, dx, match_cost, mismatch_cost, sign):
	# A permutation of moves states which string will have its character consumed (1)
	# or will have a gap instead (0).
	perms = [  # Exclude perm with all gaps. Not allowed!
		[0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 1], [0, 1, 0, 0], [0, 1, 0, 1],
		[0, 1, 1, 0], [0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 1], [1, 0, 1, 0],
		[1, 0, 1, 1], [1, 1, 0, 0], [1, 1, 0, 1], [1, 1, 1, 0], [1, 1, 1, 1]
	]
	best_score = -999999999
	best_perm = None
	for perm in perms:
		# Check if perm is applicable.
		if not all([not perm[0] or i > 0, not perm[1] or j > 0, not perm[2] or k > 0, not perm[3] or l > 0]):
			continue
		ax_p = ax if perm[0] else '-'
		bx_p = bx if perm[1] else '-'
		cx_p = cx if perm[2] else '-'
		dx_p = dx if perm[3] else '-'
		tmp = [ax_p, bx_p, cx_p, dx_p]
		score = D[i-perm[0]][j-perm[1]][k-perm[2]][l-perm[3]]  # Include corresponding previous D value.
		# Add sum-of-pairs score.
		for p0 in range(len(tmp)-1):
			for p1 in range(p0+1, len(tmp)):
				score += sign*match_cost if tmp[p0] == tmp[p1] else sign*mismatch_cost
		# Keep track of best score and corresponding permutation of moves.
		if score > best_score:
			best_score = score
			best_perm = perm
	return best_score, best_perm


match_cost = 0
mismatch_cost = 1
sign = -1
D = [[[[0 for l in range(len(d)+1)] for k in range(len(c)+1)] for j in range(len(b)+1)] for i in range(len(a)+1)]
# Remeber which permutation led to the cell at i,j,k,l.
backtrace = [[[[None for l in range(len(d)+1)] for k in range(len(c)+1)] for j in range(len(b)+1)] for i in range(len(a)+1)]
# Initialize base cases. Watch out, this only initializes 1D-lines in the multidimensional space.
# Those lines correspond to only advancing one string -> other n-1 strings get gaps -> sum-of-pairs score = n-1.
# Planes between these base case lines need to be calculated as well.
for i in range(1, len(a)+1):
	D[i][0][0][0] = i*sign*3*mismatch_cost
	backtrace[i][0][0][0] = [1, 0, 0, 0]
for j in range(1, len(b)+1):
	D[0][j][0][0] = j*sign*3*mismatch_cost
	backtrace[0][j][0][0] = [0, 1, 0, 0]
for k in range(1, len(c)+1):
	D[0][0][k][0] = k*sign*3*mismatch_cost
	backtrace[0][0][k][0] = [0, 0, 1, 0]
for l in range(1, len(d)+1):
	D[0][0][0][l] = l*sign*3*mismatch_cost
	backtrace[0][0][0][l] = [0, 0, 0, 1]
for i in range(len(a)+1):
	for j in range(len(b)+1):
		for k in range(len(c)+1):
			for l in range(len(d)+1):
				# If at one of the base case lines, skip.
				if i+j+k+l <= 1:
					continue
				best_score, best_perm = get_best_score(D, i, j, k, l, a[i-1], b[j-1], c[k-1], d[l-1],
					match_cost, mismatch_cost, sign)
				D[i][j][k][l] = best_score
				backtrace[i][j][k][l] = best_perm

i, j, k, l = len(a), len(b), len(c), len(d)
a_align = []
b_align = []
c_align = []
d_align = []
while i > 0 or j > 0 or k > 0 or l > 0:
	perm = backtrace[i][j][k][l]
	a_align.append(a[i-1] if perm[0] else '-')
	b_align.append(b[j-1] if perm[1] else '-')
	c_align.append(c[k-1] if perm[2] else '-')
	d_align.append(d[l-1] if perm[3] else '-')
	i, j, k, l = i-perm[0], j-perm[1], k-perm[2], l-perm[3]

print(D[-1][-1][-1][-1])
print(''.join(a_align[::-1]))
print(''.join(b_align[::-1]))
print(''.join(c_align[::-1]))
print(''.join(d_align[::-1]))

# Verify alignment score.
alignment_score = 0
tmp = [a_align, b_align, c_align, d_align]
for i in range(len(a_align)):
	for j in range(len(tmp)-1):
		for k in range(j+1, len(tmp)):
			if tmp[j][i] != tmp[k][i]:
				alignment_score -= 1
error_msg = "D value -> {} != {} <- score of result alignment".format(D[-1][-1][-1][-1], alignment_score)
assert D[-1][-1][-1][-1] == alignment_score, error_msg