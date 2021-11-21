import string


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

with open("rosalind_edta.txt") as file:
	dna_list = parse_FASTA(file)

s, t = [x.content for x in dna_list]
gamma = 1
match_cost = 0
mismatch_cost = 1
dists = [[d for d in range(len(s) + 1)]] + [[None for j in range(len(s) + 1)] for i in range(len(t))]
for i in range(1, len(t) + 1):
	dists[i][0] = gamma * i
	for j in range(1, len(s) + 1):
		dists[i][j] = min(dists[i - 1][j - 1] + (match_cost if s[j - 1] == t[i - 1] else mismatch_cost), dists[i - 1][j] + gamma, dists[i][j - 1] + gamma)

print(dists[-1][-1])
s_alig = []
t_alig = []
i = len(t)
j = len(s)
while i > 0 and j > 0:
	insertion_val = dists[i - 1][j] + gamma
	deletion_val = dists[i][j - 1] + gamma
	subst_val = dists[i - 1][j - 1] + (match_cost if s[j - 1] == t[i - 1] else mismatch_cost)
	min_val = min(subst_val, insertion_val, deletion_val)
	if subst_val == min_val:
		i -= 1
		j -= 1
		s_alig.append(s[j])
		t_alig.append(t[i])
	elif deletion_val == min_val:
		j -= 1
		s_alig.append(s[j])
		t_alig.append('-')
	elif insertion_val == min_val:
		i -= 1
		s_alig.append('-')
		t_alig.append(t[i])
	else:
		raise Exception()
if i > 0:
	t_alig += ['-'] * i
	i = 0
if j > 0:
	s_alig += ['-'] * j
	j = 0

print(''.join(reversed(s_alig)))
print(''.join(reversed(t_alig)))