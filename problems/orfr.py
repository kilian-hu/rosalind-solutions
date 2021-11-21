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

############
# GO HERE! #
############
with open("rosalind_orfr.txt") as file:
	dna = file.readline().strip()

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

print(max(set(find_proteins(dna) + find_proteins(dna_rc)), key=len))