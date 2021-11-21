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


with open("rosalind_gcon.txt") as file:
	dna_list = parse_FASTA(file)

a, b = [x.content for x in dna_list]

score_mat = [  # substitution penalty
	[4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2],
	[0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
	[-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3],
	[-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2],
	[-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3],
	[0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3],
	[-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2],
	[-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1],
	[-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2],
	[-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1],
	[-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1],
	[-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2],
	[-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3],
	[-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1],
	[-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2],
	[1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2],
	[0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2],
	[0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1],
	[-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2],
	[-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7]
]
score_mat_map = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9,
	'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}

# g(k) = alpha + beta * k
alpha = -5  # Gap initial penalty.
beta = 0  # Gap continuation penalty.

D = [[None for j in range(len(a)+1)] for i in range(len(b)+1)]
P = [[None for j in range(len(a)+1)] for i in range(len(b)+1)]
Q = [[None for j in range(len(a)+1)] for i in range(len(b)+1)]
D[0][0], P[0][0], Q[0][0] = 0, 0, 0
for i in range(1, len(b)+1):
	D[i][0] = alpha + i * beta
	Q[i][0] = -sys.maxsize
for j in range(1, len(a)+1):
	D[0][j] = alpha + j * beta
	P[0][j] = -sys.maxsize

for i in range(1, len(b)+1):
	for j in range(1, len(a)+1):
		P[i][j] = max(D[i-1][j] + alpha + beta, P[i-1][j] + beta)
		Q[i][j] = max(D[i][j-1] + alpha + beta, Q[i][j-1] + beta) 
		D[i][j] = max(P[i][j], Q[i][j], D[i-1][j-1] + score_mat[score_mat_map[a[j-1]]][score_mat_map[b[i-1]]])

print(D[-1][-1])