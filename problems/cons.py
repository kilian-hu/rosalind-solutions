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
nucleotides = 'ACGT'
with open('rosalind_cons.txt') as file:
	dna_list = parse_FASTA(file)

dna_list = [x.content for x in dna_list]
n = len(dna_list[0])
assert all([len(x) == n for x in dna_list])
profile = dict([(x, [0 for _ in range(n)]) for x in nucleotides])
for dna in dna_list:
	for i in range(n):
		profile[dna[i]][i] += 1

consensus = [None for _ in range(n)]
for i in range(n):
	max_val = -1
	max_val_nucleotide = None
	for x in nucleotides:
		if profile[x][i] > max_val:
			max_val = profile[x][i]
			max_val_nucleotide = x
	consensus[i] = max_val_nucleotide

print(''.join(consensus))
for x in nucleotides:
	print(x + ': ', end='')
	print (' '.join([str(c) for c in profile[x]]))