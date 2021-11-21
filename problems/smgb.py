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


############
# GO HERE! #
############
with open("rosalind_smgb.txt") as file:
	dna_list = parse_FASTA(file)

a, b = [x.content for x in dna_list]

mode = 'max'
func = {'max': max, 'min': min}[mode]
sign = {'max': -1, 'min': 1}[mode]  # When maximizing, negative is for punishment
comp_op = {'max': lambda x, y: x > y, 'min': lambda x, y: x < y}[mode] 

gamma = sign * 1  # Gap penalty.
match_cost = -sign * 1
substitution_cost = sign * 1

D = [[None for _ in range(len(a)+1)] for _ in range(len(b)+1)]
for i in range(len(b)+1):
	D[i][0] = 0
for j in range(len(a)+1):
	D[0][j] = 0

for i in range(1, len(b)+1):
	for j in range(1, len(a)+1):
		match_val = match_cost if a[j-1] == b[i-1] else substitution_cost
		a_gamma = gamma if i < len(b) else 0
		b_gamma = gamma if j < len(a) else 0
		D[i][j] = func(D[i-1][j-1]+match_val, D[i-1][j]+b_gamma, D[i][j-1]+a_gamma)

print(D[-1][-1])

a_alig = []
b_alig = []
cur_i, cur_j = len(b), len(a)
while cur_i > 0 or cur_j > 0:
	match_val = match_cost if a[cur_j-1] == b[cur_i-1] else substitution_cost
	a_gamma = gamma if cur_i not in [0, len(b)] else 0  # In forward pass, this is handled by
	b_gamma = gamma if cur_j not in [0, len(a)] else 0  # initialization. Here, check both.
	if (cur_i > 0 and cur_j > 0) and D[cur_i][cur_j] == D[cur_i-1][cur_j-1] + match_val:
		a_alig.append(a[cur_j-1])
		b_alig.append(b[cur_i-1])
		cur_i -= 1
		cur_j -= 1
	elif cur_i > 0 and D[cur_i][cur_j] == D[cur_i-1][cur_j] + b_gamma:
		a_alig.append('-')
		b_alig.append(b[cur_i-1])
		cur_i -= 1
	elif cur_j > 0 and D[cur_i][cur_j] == D[cur_i][cur_j-1] + a_gamma:
		a_alig.append(a[cur_j-1])
		b_alig.append('-')
		cur_j -= 1
	else:
		raise Exception()

print(''.join(a_alig[::-1]))
print(''.join(b_alig[::-1]))