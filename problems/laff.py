import string
import sys


VALID_CHARS = set([c for c in string.ascii_letters] + [c for c in string.digits])


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


BLOSUM62 = [  # Substitution penalty BLOSUM62
	[  4,   0,  -2,  -1,  -2,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,  -1,  -1,   1,   0,   0,  -3,  -2],
	[  0,   9,  -3,  -4,  -2,  -3,  -3,  -1,  -3,  -1,  -1,  -3,  -3,  -3,  -3,  -1,  -1,  -1,  -2,  -2],
	[ -2,  -3,   6,   2,  -3,  -1,  -1,  -3,  -1,  -4,  -3,   1,  -1,   0,  -2,   0,  -1,  -3,  -4,  -3],
	[ -1,  -4,   2,   5,  -3,  -2,   0,  -3,   1,  -3,  -2,   0,  -1,   2,   0,   0,  -1,  -2,  -3,  -2],
	[ -2,  -2,  -3,  -3,   6,  -3,  -1,   0,  -3,   0,   0,  -3,  -4,  -3,  -3,  -2,  -2,  -1,   1,   3],
	[  0,  -3,  -1,  -2,  -3,   6,  -2,  -4,  -2,  -4,  -3,   0,  -2,  -2,  -2,   0,  -2,  -3,  -2,  -3],
	[ -2,  -3,  -1,   0,  -1,  -2,   8,  -3,  -1,  -3,  -2,   1,  -2,   0,   0,  -1,  -2,  -3,  -2,   2],
	[ -1,  -1,  -3,  -3,   0,  -4,  -3,   4,  -3,   2,   1,  -3,  -3,  -3,  -3,  -2,  -1,   3,  -3,  -1],
	[ -1,  -3,  -1,   1,  -3,  -2,  -1,  -3,   5,  -2,  -1,   0,  -1,   1,   2,   0,  -1,  -2,  -3,  -2],
	[ -1,  -1,  -4,  -3,   0,  -4,  -3,   2,  -2,   4,   2,  -3,  -3,  -2,  -2,  -2,  -1,   1,  -2,  -1],
	[ -1,  -1,  -3,  -2,   0,  -3,  -2,   1,  -1,   2,   5,  -2,  -2,   0,  -1,  -1,  -1,   1,  -1,  -1],
	[ -2,  -3,   1,   0,  -3,   0,   1,  -3,   0,  -3,  -2,   6,  -2,   0,   0,   1,   0,  -3,  -4,  -2],
	[ -1,  -3,  -1,  -1,  -4,  -2,  -2,  -3,  -1,  -3,  -2,  -2,   7,  -1,  -2,  -1,  -1,  -2,  -4,  -3],
	[ -1,  -3,   0,   2,  -3,  -2,   0,  -3,   1,  -2,   0,   0,  -1,   5,   1,   0,  -1,  -2,  -2,  -1],
	[ -1,  -3,  -2,   0,  -3,  -2,   0,  -3,   2,  -2,  -1,   0,  -2,   1,   5,  -1,  -1,  -3,  -3,  -2],
	[  1,  -1,   0,   0,  -2,   0,  -1,  -2,   0,  -2,  -1,   1,  -1,   0,  -1,   4,   1,  -2,  -3,  -2],
	[  0,  -1,  -1,  -1,  -2,  -2,  -2,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,   1,   5,   0,  -2,  -2],
	[  0,  -1,  -3,  -2,  -1,  -3,  -3,   3,  -2,   1,   1,  -3,  -2,  -2,  -3,  -2,   0,   4,  -3,  -1],
	[ -3,  -2,  -4,  -3,   1,  -2,  -2,  -3,  -3,  -2,  -1,  -4,  -4,  -2,  -3,  -3,  -2,  -3,  11,   2],
	[ -2,  -2,  -3,  -2,   3,  -3,   2,  -1,  -2,  -1,  -1,  -2,  -3,  -1,  -2,  -2,  -2,  -1,   2,   7]
]
BLOSUM62_ENTRIES = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
	'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
BLOSUM62_MAP = dict(zip(BLOSUM62_ENTRIES, range(len(BLOSUM62_ENTRIES))))

PAM250 = [  # Substitution penalty PAM250
	[  2,  -2,   0,   0,  -3,   1,  -1,  -1,  -1,  -2,  -1,   0,   1,   0,  -2,   1,   1,   0,  -6,  -3],
	[ -2,  12,  -5,  -5,  -4,  -3,  -3,  -2,  -5,  -6,  -5,  -4,  -3,  -5,  -4,   0,  -2,  -2,  -8,   0],
	[  0,  -5,   4,   3,  -6,   1,   1,  -2,   0,  -4,  -3,   2,  -1,   2,  -1,   0,   0,  -2,  -7,  -4],
	[  0,  -5,   3,   4,  -5,   0,   1,  -2,   0,  -3,  -2,   1,  -1,   2,  -1,   0,   0,  -2,  -7,  -4],
	[ -3,  -4,  -6,  -5,   9,  -5,  -2,   1,  -5,   2,   0,  -3,  -5,  -5,  -4,  -3,  -3,  -1,   0,   7],
	[  1,  -3,   1,   0,  -5,   5,  -2,  -3,  -2,  -4,  -3,   0,   0,  -1,  -3,   1,   0,  -1,  -7,  -5],
	[ -1,  -3,   1,   1,  -2,  -2,   6,  -2,   0,  -2,  -2,   2,   0,   3,   2,  -1,  -1,  -2,  -3,   0],
	[ -1,  -2,  -2,  -2,   1,  -3,  -2,   5,  -2,   2,   2,  -2,  -2,  -2,  -2,  -1,   0,   4,  -5,  -1],
	[ -1,  -5,   0,   0,  -5,  -2,   0,  -2,   5,  -3,   0,   1,  -1,   1,   3,   0,   0,  -2,  -3,  -4],
	[ -2,  -6,  -4,  -3,   2,  -4,  -2,   2,  -3,   6,   4,  -3,  -3,  -2,  -3,  -3,  -2,   2,  -2,  -1],
	[ -1,  -5,  -3,  -2,   0,  -3,  -2,   2,   0,   4,   6,  -2,  -2,  -1,   0,  -2,  -1,   2,  -4,  -2],
	[  0,  -4,   2,   1,  -3,   0,   2,  -2,   1,  -3,  -2,   2,   0,   1,   0,   1,   0,  -2,  -4,  -2],
	[  1,  -3,  -1,  -1,  -5,   0,   0,  -2,  -1,  -3,  -2,   0,   6,   0,   0,   1,   0,  -1,  -6,  -5],
	[  0,  -5,   2,   2,  -5,  -1,   3,  -2,   1,  -2,  -1,   1,   0,   4,   1,  -1,  -1,  -2,  -5,  -4],
	[ -2,  -4,  -1,  -1,  -4,  -3,   2,  -2,   3,  -3,   0,   0,   0,   1,   6,   0,  -1,  -2,   2,  -4],
	[  1,   0,   0,   0,  -3,   1,  -1,  -1,   0,  -3,  -2,   1,   1,  -1,   0,   2,   1,  -1,  -2,  -3],
	[  1,  -2,   0,   0,  -3,   0,  -1,   0,   0,  -2,  -1,   0,   0,  -1,  -1,   1,   3,   0,  -5,  -3],
	[  0,  -2,  -2,  -2,  -1,  -1,  -2,   4,  -2,   2,   2,  -2,  -1,  -2,  -2,  -1,   0,   4,  -6,  -2],
	[ -6,  -8,  -7,  -7,   0,  -7,  -3,  -5,  -3,  -2,  -4,  -4,  -6,  -5,   2,  -2,  -5,  -6,  17,   0],
	[ -3,   0,  -4,  -4,   7,  -5,   0,  -1,  -4,  -1,  -2,  -2,  -5,  -4,  -4,  -3,  -3,  -2,   0,  10]
]
PAM250_ENTRIES = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
	'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
PAM250_MAP = dict(zip(PAM250_ENTRIES, range(len(PAM250_ENTRIES))))

############
# GO HERE! #
############
with open("rosalind_laff.txt") as file:
	dna_list = parse_FASTA(file)

a, b = [x.content for x in dna_list]

score_mat = BLOSUM62
score_mat_map = BLOSUM62_MAP
mode = 'max'
func = {'max': max, 'min': min}[mode]
sign = {'max': -1, 'min': 1}[mode]  # When maximizing, negative is for punishment
comp_op = {'max': lambda x, y: x > y, 'min': lambda x, y: x < y}[mode] 
# g(k) := alpha + beta * k
alpha = sign * 11  # Gap initial penalty.
beta = sign * 1  # Gap continuation penalty.

D = [[None for j in range(len(a)+1)] for i in range(len(b)+1)]
P = [[None for j in range(len(a)+1)] for i in range(len(b)+1)]
Q = [[None for j in range(len(a)+1)] for i in range(len(b)+1)]
for i in range(len(b)+1):
	D[i][0] = 0
	Q[i][0] = sign * sys.maxsize
for j in range(len(a)+1):
	D[0][j] = 0
	P[0][j] = sign * sys.maxsize
P[0][0], Q[0][0] = 0, 0

best_val = sign * sys.maxsize
best_i, best_j = None, None
for i in range(1, len(b)+1):
	for j in range(1, len(a)+1):
		P[i][j] = func(D[i-1][j] + alpha, P[i-1][j] + beta)
		Q[i][j] = func(D[i][j-1] + alpha, Q[i][j-1] + beta)
		match_score = D[i-1][j-1] + score_mat[score_mat_map[a[j-1]]][score_mat_map[b[i-1]]]
		D[i][j] = func(0, match_score, P[i][j], Q[i][j])
		if comp_op(D[i][j], best_val):
			best_val = D[i][j]
			best_i, best_j = i, j

print(best_val)

a_alig = []
b_alig = []
cur_i, cur_j = best_i, best_j
cur_mat = 'D'
while comp_op(D[cur_i][cur_j], 0):
	if cur_mat == 'D':
		match_score = D[cur_i-1][cur_j-1] + score_mat[score_mat_map[a[cur_j-1]]][score_mat_map[b[cur_i-1]]]
		if D[cur_i][cur_j] ==  match_score:
			a_alig.append(a[cur_j-1])
			b_alig.append(b[cur_i-1])
			cur_i -= 1
			cur_j -= 1
		elif D[cur_i][cur_j] == P[cur_i][cur_j]:
			cur_mat = 'P'
		elif D[cur_i][cur_j] == Q[cur_i][cur_j]:
			cur_mat = 'Q'
		else:
			raise Exception()
	elif cur_mat == 'P':
		if P[cur_i][cur_j] != P[cur_i-1][cur_j] + beta:
			cur_mat = 'D'
		b_alig.append(b[cur_i-1])
		cur_i -= 1
	elif cur_mat == 'Q':
		if Q[cur_i][cur_j] != Q[cur_i][cur_j-1] + beta:
			cur_mat = 'D'
		a_alig.append(a[cur_j-1])
		cur_j -= 1
	else:
		raise Exception()

print(''.join(a_alig[::-1]))
print(''.join(b_alig[::-1]))