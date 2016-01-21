import os
from collections import Counter


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
            self.dna = self.revcomp(self.dna)
        codons = self.splitSeqToCodons()
        amino_acids = ''
        for codon in codons:
            amino_acids = amino_acids + self.translateCodon(codon)
        return amino_acids
        
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





fasta = mydir + 'EcoGene_no_pseudo.fa'

class_test = classFASTA(fasta)
sequence_1 = class_test.readFASTA()[0][1]
print dnaTranslation(mydir,sequence_1, 1).codonFrequency()
