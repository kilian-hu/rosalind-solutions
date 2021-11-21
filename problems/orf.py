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


DNA_CODON_MAP = {
	"TTT": 'F', "TTC": 'F', "TTA": 'L',    "TTG": 'L',
	"TCT": 'S', "TCC": 'S', "TCA": 'S',    "TCG": 'S',
	"TAT": 'Y', "TAC": 'Y', "TAA": 'Stop', "TAG": 'Stop',
	"TGT": 'C', "TGC": 'C', "TGA": 'Stop', "TGG": 'W',
	"CTT": 'L', "CTC": 'L', "CTA": 'L',    "CTG": 'L',
	"CCT": 'P', "CCC": 'P', "CCA": 'P',    "CCG": 'P',
	"CAT": 'H', "CAC": 'H', "CAA": 'Q',    "CAG": 'Q',
	"CGT": 'R', "CGC": 'R', "CGA": 'R',    "CGG": 'R',
	"ATT": 'I', "ATC": 'I', "ATA": 'I',    "ATG": 'M',
	"ACT": 'T', "ACC": 'T', "ACA": 'T',    "ACG": 'T',
	"AAT": 'N', "AAC": 'N', "AAA": 'K',    "AAG": 'K',
	"AGT": 'S', "AGC": 'S', "AGA": 'R',    "AGG": 'R',
	"GTT": 'V', "GTC": 'V', "GTA": 'V',    "GTG": 'V',
	"GCT": 'A', "GCC": 'A', "GCA": 'A',    "GCG": 'A',
	"GAT": 'D', "GAC": 'D', "GAA": 'E',    "GAG": 'E',
	"GGT": 'G', "GGC": 'G', "GGA": 'G',    "GGG": 'G'
}

DNA_START_CODON = "ATG"
DNA_STOP_CODONS = [key for key, val in DNA_CODON_MAP.items() if val == "Stop"]
DNA_COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

RNA_CODON_MAP = {
	"UUU": 'F', "UUC": 'F', "UUA": 'L',    "UUG": 'L',
	"UCU": 'S', "UCC": 'S', "UCA": 'S',    "UCG": 'S',
	"UAU": 'Y', "UAC": 'Y', "UAA": 'Stop', "UAG": 'Stop',
	"UGU": 'C', "UGC": 'C', "UGA": 'Stop', "UGG": 'W',
	"CUU": 'L', "CUC": 'L', "CUA": 'L',    "CUG": 'L',
	"CCU": 'P', "CCC": 'P', "CCA": 'P',    "CCG": 'P',
	"CAU": 'H', "CAC": 'H', "CAA": 'Q',    "CAG": 'Q',
	"CGU": 'R', "CGC": 'R', "CGA": 'R',    "CGG": 'R',
	"AUU": 'I', "AUC": 'I', "AUA": 'I',    "AUG": 'M',
	"ACU": 'T', "ACC": 'T', "ACA": 'T',    "ACG": 'T',
	"AAU": 'N', "AAC": 'N', "AAA": 'K',    "AAG": 'K',
	"AGU": 'S', "AGC": 'S', "AGA": 'R',    "AGG": 'R',
	"GUU": 'V', "GUC": 'V', "GUA": 'V',    "GUG": 'V',
	"GCU": 'A', "GCC": 'A', "GCA": 'A',    "GCG": 'A',
	"GAU": 'D', "GAC": 'D', "GAA": 'E',    "GAG": 'E',
	"GGU": 'G', "GGC": 'G', "GGA": 'G',    "GGG": 'G'
}

RNA_START_CODON = "AUG"
RNA_STOP_CODONS = [key for key, val in RNA_CODON_MAP.items() if val == "Stop"]
RNA_COMPLEMENTS = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}


############
# GO HERE! #
############
with open("rosalind_orf.txt") as file:
	dna_list = parse_FASTA(file)

dna = ''.join(dna_list[0].content)
dna_rc = ''.join(DNA_COMPLEMENTS[c] for c in reversed(dna))

def find_proteins(s):
	starts = []
	for i in range(len(s)-2):
		if s[i:i+3] == DNA_START_CODON:
			starts.append(i)
	proteins = []
	for start in starts:
		cur_protein = []
		for i in range(start, len(s)-(len(s)-start)%3, 3):
			aa = DNA_CODON_MAP[s[i:i+3]]
			if aa == "Stop":
				proteins.append(''.join(cur_protein))
				break
			cur_protein.append(aa)
	return proteins

for p in set(find_proteins(dna) + find_proteins(dna_rc)):
	print(p)