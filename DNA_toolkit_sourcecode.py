import collections
from collections import Counter
# from colorama import init, Fore; init(autoreset=True)
import colorama; from colorama import Style, Fore, Back, init; colorama.init(autoreset=True)
import random


##############################################################################################################################################################
#####   Below are some constants used in the functions   #####################################################################################################
##############################################################################################################################################################

Nucleotides     = ['A', 'C', 'G', 'T']
DNA_Nucleotides = ['A', 'C', 'G', 'T']

DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
Transcription_dict = {'G': 'C', 'C': 'G', 'T': 'A', 'A': 'U'}
DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

RNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUA": "I", "AUU": "I", "AUC": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "UAA": "_", "UAG": "_", "UGA": "_"
}

##############################################################################################################################################################
#####   Now functions themselves   ###########################################################################################################################
##############################################################################################################################################################

def validateSeq(dna_seq):
    """
    Check string to check that it's a DNA string, i.e. it only contains 'A', 'T', 'C', 'G', 
    either capital or lowercase.
    """
    tmpseq = dna_seq.upper()
    temp = 0
    for nuc in tmpseq:
        if nuc in Nucleotides:
            temp += 0
        if nuc not in Nucleotides:
            temp += 1
    if temp == 0:
        #print("Sequence consists only of nucleotides")
        #print(tmpseq)
        return True
    if temp != 0:
        #print("FALSE: not a nucleotide sequence")
        return False

# check that strings over multiple lines in a file 'input.txt' are DNA sequences
# print only those strings that are DNA sequences into output file
def validateSequences_fromFile(a, b):
    """
    Check that lines in input file 'input.txt' are DNA sequences; \n
    print only those strings that are DNA sequences into output file; \n
    ```validateSequences_fromFile('input.txt', 'output.txt')```
    """
    with open(a, 'r') as input, open(b, 'w') as output:
        temp1 = input.readlines()
        temp2 = [item.replace("\n", "") for item in temp1]
        temp3 = [item.replace(" ", "") for item in temp2]
        bin = 0
        for line in temp3:
            if validateSeq(str(line)) == True:
                output.write(line.upper() + "\n")
            if validateSeq(str(line)) == False:
                print('WARNING! incorrect sequence encountered. The sequence will be skipped')
                bin += 1
                pass
        print(f"\nProgram successful! \n{bin} sequences have been skipped")

def transcription(sequence):
    """returns reverse complement (RNA) of DNA sequence"""
    return "".join([Transcription_dict[nuc] for nuc in sequence])[::-1]

def translation_0(seq, init_pos=0):
    """Translates a DNA sequence into an aminoacid sequence"""
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]

def translation_1(sequence, k=3):
    """Translates RNA sequence into AA sequence"""
    res = ""
    for i in range(0, len(sequence) -k +1, k):
        subseq = sequence[i: i+k]
        if DNA_Codons[subseq] != "_":
            res += DNA_Codons[subseq]
    return res

def ds_dna_maker(sequence):
    """Outputs sequence along with its complementary strand"""
    # print(f"DNA String + Reverse Complement:\n5' {sequence} 3'")
    print(f"5' {sequence} 3'")
    print(f"   {''.join(['|' for c in range(len(sequence))])}")
    print(f"3' {reverse_complement_2(sequence)[::-1]} 5'")    

def T2U_replace(seq):
    """Replace all 'T' with 'U' nucleotides"""
    return seq.replace('T', 'U')

def codon_usage(seq, aminoacid):
    """returns counts of specified amino acid in the DNA sequence"""
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i+3]] == aminoacid:
            tmpList.append(seq[i:i+3])
    #print(tmpList)
    freqDict = {codon: tmpList.count(codon) for codon in tmpList}
    return freqDict
#print( codon_usage("ATAATTATCTGAATAATGATATTAATGAATCCTCGTTAA", "I"))

def gen_reading_frames(seq):
    """Generate the six reading frames of a DNA sequence, 
    including the reverse complement"""
    frames = []
    frames.append(translation_0(seq, 0))
    frames.append(translation_0(seq, 1))
    frames.append(translation_0(seq, 2))
    frames.append(translation_0(reverse_complement_0(seq), 0))
    frames.append(translation_0(reverse_complement_0(seq), 1))
    frames.append(translation_0(reverse_complement_0(seq), 2))
    return frames
#print(gen_reading_frames(seq))
#for frame in gen_reading_frames("ATGGACATGCAGTAGCAGTAGCATCAA"):
#    print(frame)

def rabbits_fib(n):
    """ How many pairs of rabbits are present after n months, given that each M+F pair produces exactly 2 offspring; 
    from Rosalind (http://rosalind.info/problems/fib/), k = 2 """
    if n == 1 or n == 2:
        return 1
    a = [None] * (n+1)
    b = [None] * (n+1)
    a[1] = 1
    a[2] = 0
    b[1] = 0
    b[2] = 1
    for i in range(3, n+1):
        a[i] = a[i-1] + a[i-2]
        b[i] = b[i-1] + b[i-2]
    return a[n] + b[n]
# print( rabbits_fib(7) )

def rabbits_nk(n, k):
    """shows how many pairs of rabbits there are after n months;
    problem from Rosalind 'Rabbits and Recurrence Relations' 
    (http://rosalind.info/problems/fib/), k = any that you choose"""
    if n == 1 or n == 2:
        return 1
    b = [None] * (n+1)
    m = [None] * (n+1)
    b[1] = 1
    b[2] = 0
    m[1] = 0
    m[2] = 1
    for i in range(3, n+1):
        b[i] = m[i-1] * k
        m[i] = b[i-1] + m[i-1]
    return b[n] + m[n]
#print (rabbits_nk(33, 5) )

def proteins_from_rf(aa_seq):
    """Compute all possible proteins in an AA seq and return a list of possible proteins"""
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            # STOP accumulating amino acids if _ - STOP was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            # START accumulating amino acids if M - START was found
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins
# test_rf_frame = ['L', 'M', 'A', 'M', 'T', 'A', 'L', 'V', 'V', 'L', 'V', 'R', 'R', 'G', 'S', 'V', 'G', '_', 'H']
# print( proteins_from_rf(test_rf_frame) )

def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    """Compute all possible proteins for all open reading frames"""
    """Protine Search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
    """API can be used to pull protein info"""
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startReadPos: endReadPos])
    else:
        rfs = gen_reading_frames(seq)
    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)
    if ordered:
        return sorted(res, key=len, reverse=True) # sort by length; reverse = from longest to shortest; 
    return res

def find_motif(s, t):
	"""Finds at which positions (if string starts at 1) of string 's' you can find substring 't'
    """
	i = 0
	listname = []
	while i <= len(s) - len(t):
		if t == s[i:i+len(t)]:
			listname.append(i+1)
		i = i + 1
	return listname






# Print reverse complement from ONE input sequence. Uses translation table - pythonic approach. A little bit faster. 
# ver1
def reverse_complement_0(seq):
    mapping = str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)[::-1]

# Print reverse complement from ONE input sequence
# ver2
def reverse_complement_1(sequence):
    print(sequence)
    rev_seq = [i for i in sequence[::-1]]
    print(f"rev_seq: {rev_seq}")
    rev_seq2 = [item.upper() for item in rev_seq]
    print(f"rev_seq2: {rev_seq2}")
    temp = []
    for item in rev_seq2:
        if item == "C":
            temp.append("G")
        if item == "G":
            temp.append("C")
        if item == "A":
            temp.append("T")
        if item == "T":
            temp.append("A")
    print(f"temp: {temp}")
    print("".join(temp))

# Print reverse complement from ONE input sequence
# ver3
def reverse_complement_2(sequence):
    """this is an example of a docstring - explanation of the function boi"""
    return "".join([DNA_ReverseComplement[nuc] for nuc in sequence])[::-1]



# Read lines of sequences from 'input.txt', output reverse complements into 'output.txt'
def reverse_complement_3(a, b):
    """
    Flips the sequence and substitutes nucleotides
    with their complementary counterparts. 
    """
    validateSequences_fromFile(a, b)
    with open(b, 'r') as input:
        temp1 = input.readlines()
        temp2 = [item.replace("\n", "") for item in temp1]
        print(f"temp2: {temp2}")
        with open(b, 'w') as output:
            for sequence in temp2:
                print(sequence)
                rev_seq = [i for i in sequence[::-1]]
                print(f"rev_seq: {rev_seq}")

                for item in rev_seq:
                    if item == "C":
                        output.write("G")
                    if item == "G":
                        output.write("C")
                    if item == "A":
                        output.write("T")
                    if item == "T":
                        output.write("A")
                output.write("\n")


def max_gc_content_0(filename):
	"""Reads file with sequences in .fasta format,
	returns the highest %GC value from the sequences\n
    Example: max_gc_content_0('rosalind_gc.txt')"""
    
	# Read FASTA file and store lines as elements of a list
	def readFile(filePath):
		"""Read file, return a list of lines"""
		with open(filePath, 'r') as f:
			return [l.strip() for l in f.readlines()]
	FASTAFile = readFile('rosalind_gc.txt')
	# Dictionary for labels + data
	FASTADict = {}
	# string for holding the current label
	FASTALabel = ""
	# clean and prepare our data
	# convert FASTA/list file data into dictionary:
	for line in FASTAFile:
		if '>' in line:
			FASTALabel = line
			FASTADict[FASTALabel] = ""
		else:
			FASTADict[FASTALabel] += line
	# using dictionary comprehension to generate a new dictionary with GC content
	# format data: store the data in a convenient way + store GC contents
	RESULTDict = {key: gc_content_1(value) for key, value in FASTADict.items()}
	# final step - looking for max value in values of dictionary
	MaxGCKey = max(RESULTDict, key=RESULTDict.get)
	#print(f"{MaxGCKey[1:]}\n{RESULTDict[MaxGCKey]}")
	return RESULTDict[MaxGCKey]

def max_gc_content_1(filename):
	"""Reads file with sequences in .fasta format,
	returns the highest %GC value from the sequences\n
    Example: max_gc_content_1('rosalind_gc.txt')"""

	# Converts each line from input file into strings (elements) of a list 
	with open('rosalind_gc.txt', 'r') as input:
		lines = [line.strip() for line in input.readlines()]

	# Convert list of strings into a dictionary
	# e.g. {'>Key1': 'CGCGATG', '>Key2': 'AGACGAT'}
	dict1 = {}
	keyname = ""

	for element in lines:
		if '>' in element:
			keyname = element
			dict1[keyname] = ""
		else:
			dict1[keyname] = dict1[keyname] + element


	# Create a dictionary with %GC counts for each key
	dict2 = {}
	for key, value in dict1.items():
		dict2[key] = gc_content_1(value)
	
	max_gc_key = max(dict2, key=dict2.get)

	#print(f"{max_gc_key[1:]}\n{dict2[max_gc_key]}")

	return dict2[max_gc_key]



# tutorial version; returns e.g. 10
def gc_content_1(seq):
    """GC content in a DNA/RNA sequence"""
    return round( ( seq.count('C') + seq.count('G'))*100 / len(seq) )

def gc_content_11(seq):
	return float(seq.count("G") + seq.count("C"))/len(seq) * 100

# My version; returns e.g. 10.00
def gc_content_2(sequence):
    """Returns %GC content in a sequence"""
    g_content = sequence.count("G")
    c_content = sequence.count("C")
    gc_content = g_content + c_content
    #print( f"Amount of G & C nucleotides: {gc_content}" )
    percent_gc = round( gc_content * 100 / len(sequence) , 6)
    #print( f"% GC: {percent_gc}%")
    return percent_gc

def gc_content_3(seq):
	gc = [B for B in seq.upper() if B in 'GC']
	return float(len(gc))/len(seq) * 100



def gc_content_subseq(seq, k=5):
    """GC content in a DNA/RNA sub-sequence length, k (by default k=20)"""
    res = []
    for i in range(0, len(seq) -k +1, k):
        subseq = seq[i: i+k]
        res.append(gc_content_3(subseq))
    return res




# my version
def hamming_distance_calculator_1(filename):
	"""Given an input from .txt file, where there are two equal-length sequences on the first and second line, 
	calculate and return the Hamming distance between them.\n
	Hamming distance: the minimum number of symbol substitutions required to change one string into another of equal length.\n
	Example: hamming_distance_calculator('rosalind_hamm.txt')"""
	with open(filename, 'r') as input:
		lines = [line.strip() for line in input.readlines()]
		
		seq1 = lines[0]
		seq2 = lines[1]

		seq1List = [char for char in seq1]
		seq2List = [char for char in seq2]

		count = 0
		i = 0
		while i < len(seq1):
			if seq1List[i] != seq2List[i]:
				count += 1
			i += 1
		
		return count



# youtube tutorial - loop approach
def hamming_distance_calculator_2(str_1, str_2):
    h_distance = 0
    for position in range(len(str_1)):
        if str_1[position] != str_2[position]:
            h_distance += 1
    return h_distance


# youtube tutorial - using sets
def hamming_distance_calculator_3(str_1, str_2):
    nucleotide_set_1 = set([(x, y) for x, y in enumerate(str_1)])
    nucleotide_set_2 = set([(x, y) for x, y in enumerate(str_2)])
    
    # check the stored values by printing with for loop:
    print(nucleotide_set_1)
    for x in range(len(nucleotide_set_1)):
        print(sorted(nucleotide_set_1)[x], sorted(nucleotide_set_2)[x])

    return len(nucleotide_set_1.difference(nucleotide_set_2))


# youtube tutorial - using zip function
def hamming_distance_calculator_4(str_1, str_2):
    zipped_dna = zip(str_1, str_2)
    
    return  len( [(n1, n2) for n1,n2 in zipped_dna if n1 != n2] )
    
    # can do in one line:
    # return len([(n1, n2) for n1, n2 in zip(str_1, str_2) if n1 != n2])






# Count nucleotide frequencies in a sequence
# ver1
def countNucFrequency(seq):
    tmpFreqDict = {'A':0, 'C':0, 'G':0, 'T':0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

# Count nucleotide frequencies in a sequence
# ver2
def countNucFrequency_2(seq):
    return dict(collections.Counter(seq))








# my version
def codon_usage(seq, aminoacid):
    """returns counts of specified amino acid in the DNA sequence"""
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i+3]] == aminoacid:
            tmpList.append(seq[i:i+3])
    #print(tmpList)

    
    freqDict = {codon: tmpList.count(codon) for codon in tmpList}
    return freqDict

    #freqDict = dict(Counter(tmpList))
    #print(freqDict)

    #totalWight = sum(freqDict.values())
    #for seq in freqDict:
    #    freqDict[seq] = round(freqDict[seq] / totalWight, 2)
    #return freqDict
#print( codon_usage("ATAATTATCTGAATAATGATATTAATGAATCCTCGTTAA", "I"))


def gen_reading_frames(seq):
    """Generate the six reading frames of a DNA sequence, 
    including the reverse complement"""
    frames = []
    frames.append(translation_0(seq, 0))
    frames.append(translation_0(seq, 1))
    frames.append(translation_0(seq, 2))
    frames.append(translation_0(reverse_complement_0(seq), 0))
    frames.append(translation_0(reverse_complement_0(seq), 1))
    frames.append(translation_0(reverse_complement_0(seq), 2))
    return frames
#print(gen_reading_frames(seq))
#for frame in gen_reading_frames("ATGGACATGCAGTAGCAGTAGCATCAA"):
#    print(frame)



def rabbits_fib(n):
    """ How many pairs of rabbits are present after n months, given that each M+F pair produces exactly 2 offspring; 
    from Rosalind (http://rosalind.info/problems/fib/), k = 2 """
    if n == 1 or n == 2:
        return 1
    a = [None] * (n+1)
    b = [None] * (n+1)
    a[1] = 1
    a[2] = 0
    b[1] = 0
    b[2] = 1
    for i in range(3, n+1):
        a[i] = a[i-1] + a[i-2]
        b[i] = b[i-1] + b[i-2]
    return a[n] + b[n]

# print( rabbits_fib(1) )
# print( rabbits_fib(2) )
# print( rabbits_fib(3) )
# print( rabbits_fib(4) )
# print( rabbits_fib(5) )
# print( rabbits_fib(6) )
# print( rabbits_fib(7) )



def rabbits_nk(n, k):
    """shows how many pairs of rabbits there are after n months;
    problem from Rosalind 'Rabbits and Recurrence Relations' 
    (http://rosalind.info/problems/fib/), k = any that you choose"""
    if n == 1 or n == 2:
        return 1
    b = [None] * (n+1)
    m = [None] * (n+1)
    b[1] = 1
    b[2] = 0
    m[1] = 0
    m[2] = 1
    for i in range(3, n+1):
        b[i] = m[i-1] * k
        m[i] = b[i-1] + m[i-1]
    return b[n] + m[n]
#print (rabbits_nk(33, 5) )


def proteins_from_rf(aa_seq):
    """Compute all possible proteins in an AA seq and return a list of possible proteins"""
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            # STOP accumulating amino acids if _ - STOP was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            # START accumulating amino acids if M - START was found
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins


# test_rf_frame = ['L', 'M', 'A', 'M', 'T', 'A', 'L', 'V', 'V', 
# 				 'L', 'V', 'R', 'R', 'G', 'S', 'V', 'G', '_', 'H']
# print(test_rf_frame)

# print( proteins_from_rf(test_rf_frame) )



def find_motif(s, t):
	"""Finds at which positions (if string starts at 1) of string 's' you can find substring 't'"""
	i = 0
	listname = []

	while i <= len(s) - len(t):
		if t == s[i:i+len(t)]:
			listname.append(i+1)
		i = i + 1

	return listname



##############################################################################################################################################################
###   Case study of using the functions   ####################################################################################################################
##############################################################################################################################################################


if __name__ == '__main__':
    DATE = "16.august.2021"
    NAME = "A collection of different functions related to bioinformatics."
    # Create a random DNA sequence for testing: 
    # randDNAStr = ''.join([random.choice(Nucleotides) for nuc in range(50)]); print(randDNAStr)
    dnaseq = "CCGATGTTGAACTCGATACGGCTTCATTCCTACTTAACTTAGCAACATTCGCGCTA"
    print('-' * 70)
    print(f"[1] {Fore.GREEN}{Style.BRIGHT}Valid DNA sequence: {validateSeq(dnaseq)}")
    print(f"[2] {Fore.GREEN}{Style.BRIGHT}Nucleotide frequencies:"); 
    print(f"{countNucFrequency(dnaseq)}")
    print(f"[3] {Fore.GREEN}{Style.BRIGHT}DNA string and its reverse complement:"); print(f"{ds_dna_maker(dnaseq)}")
    print(f"[4] {Fore.GREEN}{Style.BRIGHT}Reverse complement:")
    print(f"{reverse_complement_0(dnaseq)}")
    print(f"[5] {Fore.GREEN}{Style.BRIGHT}GC content: {gc_content_1(dnaseq)}")
    print(f"[6] {Fore.GREEN}{Style.BRIGHT}Translation:") 
    print(f"{translation_1(dnaseq)}")
    print(f"[8] {Fore.GREEN}{Style.BRIGHT}Codon frequency (L):") 
    print(f"{codon_usage(dnaseq, 'L')}")
    print(f"[9] {Fore.GREEN}{Style.BRIGHT}Reading frames: ")
    i = 1
    for frame in gen_reading_frames(dnaseq):
        print("    " + str(i) +". " + ''.join(frame))
        i += 1
    print(f"[10] {Fore.GREEN}{Style.BRIGHT}All proteins in 6 ORFs: ")
    for prot in all_proteins_from_orfs(dnaseq, 0, 0, True):
        print(prot)
    print('-'*70)

