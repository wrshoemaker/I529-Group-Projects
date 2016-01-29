import os, argparse, random, math
from collections import Counter
from itertools import product


mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'

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


# converted into codon usage table format
def CodonTableFormat(codon_usages, **kw):
    codon_num = 0
    out = ''
    nucleotides = ['A','T','C','G']
    key = [''.join(i) for i in product(nucleotides, repeat = 3)]
    for i in range(len(key)):
        if key[i] not in codon_usages.keys():
            codon_usages[key[i]] = 0.0
    for codon in sorted(codon_usages):
        freq = "%.2f" % (round(codon_usages[codon]*100,2)) + "%"
        if codon_num == 0:
            out += codon+": "+freq
        elif codon_num % 4 != 0:
            out += "\t"+codon+": "+freq
        else:
            out += "\n"+codon+": "+freq
        codon_num += 1
    return out

# simulate randomized sequence
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

#compare the probabilty to produce the winow sequence by traning or random model
def WindowExtract(seq,size,train,random, reading_frame):
    output=''
    for i in range(len(seq)-size+1):
        window=seq[i:i+size]
        trainp=1
        randomp=1
        for j in range(int(size)/3):
            codon = window[3*j:3*j+3]
            trainp *= train[codon]
            randomp *= random[codon]
        if trainp > randomp:
            output+=str(i+1)+'\t'+str(i+size)+'\t'+str(trainp)+'\n'
    return output



# probability model with WindowSize = 99bps
def LikelihoodMode(myseq, frame, codon_usages, random_usages, threshold, RevComp = False, **kw):
    WindowSize = 99
    if RevComp == True:
        myseq = revcomp_quick_fix(myseq)
    for i in range(frame-1,len(myseq)-WindowSize, 3):
        Pc, Po, Ratio = 1.0, 1.0, 0.0
        for j in range(i, i+WindowSize, 3):
            curr_codon = myseq[j]+myseq[j+1]+myseq[j+2]
            Pc *= codon_usages[curr_codon]
            Po *= random_usages[curr_codon]
        Ratio = math.log(Pc/Po)
        #print Ratio
        if Ratio > threshold:  # return relative likelihood when it is larger than threshold
            yield str(i)+"\t"+str(i+WindowSize-1)+"\t"+str(Ratio)



#### ====== read input FASTA file and write output files ======
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--fasta_file', required=True)
    parser.add_argument('-o', '--out_table', required=True)
    parser.add_argument('-c', '--codon_table')
    parser.add_argument('-r', '--random_codon_table')
    parser.add_argument('-t','--threshold_likelihood', default=10)
    parser.add_argument('-f','--reading_frame', type = int, default=1)
    #parser.add_argument('-rc','--reverse_compliment', type = bool, default=False)
    parser.add_argument('-rc', dest='feature', action='store_true')
    parser.set_defaults(feature=False)
    args = parser.parse_args()
    fasta = mydir + 'EcoGene_no_pseudo.fa'

    class_test = classFASTA(fasta)
    ######### This command returns the nested list containing sequence names
    ######### and sequences to a flat list containing only sequences
    sequences = [ x[1] for x in class_test.readFASTA() ]
    ######### concatenated the sequences
    concatenated_seq = ''.join(sequences)
    args.feature = bool(args.feature)
    # tranform the codon usage dictionary to an optional output table
    codon_usages = dnaTranslation(mydir,concatenated_seq, args.reading_frame, reverseTL = args.feature).codonFrequency()
    codon_table = CodonTableFormat(codon_usages)

    # simulate a randomized sequence based on known nucleotide frequency
    myseq = dnaTranslation(mydir,concatenated_seq, args.reading_frame, reverseTL = args.feature).dna
    random_seq = SeqSimulator(myseq)
    random_usages = dnaTranslation(mydir,random_seq, args.reading_frame, reverseTL = args.feature).codonFrequency()
    random_table = CodonTableFormat(random_usages)

    ## generate codon usage table if required
    if 	args.codon_table:
        CodonTable = open(args.codon_table, 'w')
        CodonTable.write(codon_table)

    ## generate randomized codon table if required
    if 	args.random_codon_table:
        RandomTable = open(args.random_codon_table, 'w')
        RandomTable.write(random_table)


    #print(codon_usages['ATG'])
    ## generate the extract window
    #extract_window=WindowExtract(concatenated_seq,99,codon_usages,random_usages)
    #print(extract_window)

    # decide whether we need the negative likelihood, default is to ignore
    threshold = float(args.threshold_likelihood)

    ## read the sample fasta file
    infile = classFASTA(args.fasta_file)
    testSeq = infile.readFASTA()[0][1]

    result = LikelihoodMode(testSeq, args.reading_frame, codon_usages, random_usages, threshold, RevComp = args.feature)
    outfile = open(args.out_table, "w")
    outfile.write("start\tend\tlog(Pc/Po)\n")
    for key in result:
    	outfile.write(key+"\n")

    ## for this sample sequence, I know there is one gene from 578-992
