import os, argparse, random, math
from collections import Counter
from itertools import product


mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir) + '/'

class classFASTA:

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print "Not in FASTA format."

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list

class dnaTranslation:

    '''This class takes the fasta in as a nested list object. \
        This object is created using the above class classFASTA. \
        Codon tables are imported and are in the location  \
        ../data/codon.txt \
        AA sequence is returned as a string.'''

    def __init__(self, myDir, dna, frame, reverseTL = False):
        self.myDir = myDir + 'codon.txt'
        self.frame = frame
        self.reverseTL = reverseTL
        self.dna = dna

    def importAAtable(self):
        table = {}
        aa_list = [[item.strip() for item in line.rstrip('\r\n').split('= ', 1)[-1]] \
            for line in open(self.myDir)]
        aa_list_zip = map(list, zip(*aa_list))
        for x in aa_list_zip:
             x[2:5] = [''.join(x[2:5])]
             table[str(x[2])] = str(x[0])
        return table

    def translateCodon(self, codon):
        AAtable = self.importAAtable()
        return AAtable.get(codon.upper(), 'x')

    def splitSeqToCodons(self):
        codons = []
        for i in range(self.frame - 1, len(self.dna)-2, 3):
            codon = self.dna[i:i+3]
            codons.append(codon)
        return codons

    def revcomp(self):
        bases = 'ATGCTACG'
        comp_dict = {bases[i]:bases[i+4] for i in range(4)}
        seq = reversed(self.dna)
        list_rev_comp = [comp_dict[base] for base in seq]
        return ''.join(list_rev_comp)

    def translateDnaFrame(self):

        ''''Translates a dna sequence of a specified frame'''
        if self.reverseTL == True:
            self.dna = self.revcomp()
        codons = self.splitSeqToCodons()
        return codons
        #amino_acids = ''
        #for codon in codons:
        #    amino_acids = amino_acids + self.translateCodon(codon)
        #return amino_acids

    #def translateDnaAllFrames(self, dna):
    #    '''Translates dna sequence in 3 forward frames. \
    #        Translates from the opposite end if reverseTL == True \
    #        Returns a list of of length 6 containing the translated dna. \
    #        The first 3 items are the translated DNA sequences from frames \
    #        1, 2, and 3, respectively. Items 4, 5, and 6 contain the reverse \
    #        compliment translations for frames 1, 2, and 3, respectively'''
    #    all_translations = []
    #
    #    for frame in range(1,4):
    #        all_translations.append(self.translateDnaFrame(dna))
    #    if self.reverseTL:
    #        dnaRC = self.revcomp(dna)
    #        for frame in range(1,4):
    #            all_translations.append(self.translateDnaFrame(dnaRC))
    #    return all_translations
    #
    def codonFrequency(self):
        dnaCodons = self.translateDnaFrame()
        codonCount = Counter(dnaCodons)
        codonTotal = sum(codonCount.values(), 0.0)
        for key in codonCount:
            codonCount[key] /= codonTotal
        return codonCount


## converted into codon usage table format
## converted into codon usage table format
## converted into codon usage table format

## simulate randomized sequence
def SeqSimulator(myseq):
    randomized_seq = list(myseq)
    random.shuffle(randomized_seq)
    Simulated_sequence = ''.join(randomized_seq)
    return Simulated_sequence

def revcomp_quick_fix(dna):
    bases = 'ATGCTACG'
    comp_dict = {bases[i]:bases[i+4] for i in range(4)}
    seq = reversed(dna)
    list_rev_comp = [comp_dict[base] for base in seq]
    return ''.join(list_rev_comp)
	
def start_codon_freq(myseq_list, **kw):
    nucleotides = ['A','T','C','G']
    codons = [''.join(i) for i in product(nucleotides, repeat = 3)]
    seq_num = len(myseq_list)
    start_codons = {}
    for key in codons: 
        start_codons[key] = 0
    for key in myseq_list:
        start = key[:3]
        start_codons[start] += 1
    for key in start_codons:
        if start_codons[key] == 0:
            start_codons[key] += 1
            seq_num += 1
    for key in start_codons:
        start_codons[key] = float(start_codons[key])/seq_num
    return start_codons

def Markov_chain(myseq, **kw):
    nucleotides = ['A','T','C','G']	
    codons = [''.join(i) for i in product(nucleotides, repeat = 3)]
    codon_trans = ['-'.join(i) for i in product(codons, repeat = 2)]
    markov_table = {}
    trans_num = 0
    for key in codon_trans: 
        markov_table[key] = 0
    for x in range(3,len(myseq)-3,3):
        curr = str(myseq[x-3:x])+'-'+str(myseq[x:x+3])
        trans_num += 1
        markov_table[curr] += 1
    for key in markov_table:
        if markov_table[key] == 0:
            markov_table[key] += 1
            trans_num += 1
    for key in markov_table:
        markov_table[key] = float(markov_table[key])/trans_num
    return markov_table	

## probability model with WindowSize = 99bps
def MarkovMode(myseq, real_start, random_start, real_Markov, random_Markov, **kw):
    for x in range(0,len(myseq)-3, 3):
        codon = myseq[x:x+3]
        if x < 3:
            Pc = real_start[codon]
            Po = random_start[codon]
            Ratio = Pc/Po
            continue
        prior = myseq[x-3:x]
        trans = prior+'-'+codon
        if 'N' in trans: continue
        Pc = real_Markov[trans]
        Po = random_Markov[trans]
        Ratio *= (Pc/Po)
    Ratio = math.log(Ratio)
    return Ratio	


	
#### ====== read input FASTA file and write output files ======
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--fasta_file', required=True)
    parser.add_argument('-o', '--out_table', required=True)
    parser.add_argument('-c', '--codon_table')
    parser.add_argument('-r', '--random_codon_table')
    parser.add_argument('-t','--threshold_likelihood', default=3)
    parser.add_argument('-f','--reading_frame', type = int, default=1)
    parser.add_argument('-rc', '--reverse_compliment', dest='feature', action='store_true')
    parser.set_defaults(feature=False)
    args = parser.parse_args()
    fasta = mydir + 'ecoli-genes.fasta'

    class_test = classFASTA(fasta)
    ######### This command returns the nested list containing sequence names
    ######### and sequences to a flat list containing only sequences
    sequences = [ x[1] for x in class_test.readFASTA() ]
    sequences = [x.upper() for x in sequences]
    ######### concatenated the sequences
    concatenated_seq = ''.join(sequences[:1000])
    
    startCodons = start_codon_freq(sequences)
    Markov_table = Markov_chain(concatenated_seq)
	
    # simulate a randomized sequence based on known nucleotide frequency
    myseq = dnaTranslation(mydir,concatenated_seq, args.reading_frame, reverseTL = args.feature).dna
    random_seqs = [ SeqSimulator(x) for x in sequences ]
    random_concatseq = ''.join(random_seqs)
    randomStart = start_codon_freq(random_seqs)
    randomMarkov = Markov_chain(random_concatseq)

    ## read the sample fasta file
    infile = classFASTA(args.fasta_file)
    outfile = open(args.out_table, "w")
    outfile.write('Contig\tLog(Pc/Po) \n') 
	
    all, positive = 0, 0 	
    cutoff = args.threshold_likelihood
	
    for x in range(len(infile.readFASTA())):	
		testSeq = infile.readFASTA()[x][1]
		testID = (infile.readFASTA()[x][0]).split()[0]
		testSeq = testSeq.upper()	
		all += 1		
		results = MarkovMode(testSeq, startCodons, randomStart, Markov_table, randomMarkov)
		if results > cutoff:
			positive += 1
		outfile.write(testID+'\t'+str(results)+'\n')   

    outfile.close()
    with file(args.out_table, 'r') as original: data = original.read()
    with file(args.out_table, 'w') as modified: modified.write("Number of investigated contigs: "+str(all)+"\nNumber of contigs with log(Pc/Po)>"+str(cutoff)+": "+str(positive)+"\n\n" + data)





