from __future__ import division
import pandas as pd
from scipy.stats import gaussian_kde
import numpy as np
import  matplotlib.pyplot as plt
from collections import Counter


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


IN = "../data/test.txt"
OUT = "../data/outfile.txt"

IN_class = classFASTA(IN)

IN_seqs = [ x[1].upper() for x in IN_class.readFASTA() ]
IN_labels  = [ x[0] for x in IN_class.readFASTA() ]

IN_seqs = IN_seqs[:-1]
IN_location = [x[int(len(x)/2):] for x in IN_seqs]

OUT_class = classFASTA(OUT)

OUT_seqs = [ x[1].upper() for x in OUT_class.readFASTA() ]
OUT_labels  = [ x[0] for x in OUT_class.readFASTA() ]
OUT_last = OUT_seqs[-1]
OUT_last = OUT_last[int(len(OUT_last)/2):]
OUT_last_zipped = zip(OUT_last, OUT_last)
OUT_seqs = OUT_seqs[:-1]
OUT_location = [x[int(len(x)/2):] for x in OUT_seqs]
#print OUT_location
#print OUT_last
#OUT_location = OUT_location.append(OUT_last)
#IN_location = IN_location.append(OUT_last)

TP_not_membrane = len(OUT_last)

out_str = ''.join(OUT_location)
in_str = ''.join(IN_location)

zipped_in_out =  zip(in_str,out_str)
zipped_in_out = zipped_in_out + OUT_last_zipped
counter = 0
N = len(zipped_in_out)


counting =  Counter(zipped_in_out)

for key, value in counting.items():
    counting[key] = value / N

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

plt.bar(range(len(counting)), counting.values(), align='center')
plt.xticks(range(len(counting)), counting.keys(), fontsize = 8)


plt.tight_layout()
output = "../figures/output_plot.png"
plt.savefig(output, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.xscale()
plt.close()
