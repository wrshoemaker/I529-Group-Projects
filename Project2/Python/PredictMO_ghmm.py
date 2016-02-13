#parse the given training set protein sequence file and generate two list of protein sequence and corresponding feature of i,m,o

import sys
from itertools import groupby

#parse the training protein sequence and features
def parse_training(training_file):
	seq_set=[]
	feature_set=[]
	for line in training_file:
		tab = line.strip().split('\t')
		seq_set.append(tab[0])
		feature_set.append(tab[1])
	return seq_set,feature_set

#table of peptide length frequency in both membrane and non-membrane domain
def generate_length(feature_set):
	m_length={}
	i_length={}
	o_length={}
	for sequence in feature_set:
		tmp=[''.join(g) for k,g in groupby(sequence)]
		for each in tmp:
			if 'm' in each.lower():
				if len(each) not in m_length.keys():
					m_length[len(each)] = 1
				else:
					m_length[len(each)] += 1
			elif 'i' in each.lower():
				if len(each) not in i_length.keys():
					i_length[len(each)] = 1
				else:
					i_length[len(each)] += 1
			elif 'o' in each.lower():
				if len(each) not in o_length.keys():
					o_length[len(each)] = 1
				else:
					o_length[len(each)] += 1

	return m_length,i_length,o_length

#calculate the frequency of length for i,m,o domain
def lengthFrequency(length_set):
	lengthFreq={}
	domain_total=sum(length_set.values(),0.0)
	for key in length_set:
		lengthFreq[key]= length_set[key]/domain_total
	return lengthFreq

###======read traing file and generate GMMM model ======
if __name__ == "__main__":
	with open("../data/TMseq.ffa","r") as training_file:
		seq_set,feature_set = parse_training(training_file)
		m_length,i_length,o_length = generate_length(feature_set)
		m_freq = lengthFrequency(m_length)
		i_freq = lengthFrequency(i_length)
		o_freq = lengthFrequency(o_length)

	with open("../data/mem_length.txt","w") as mfile:
		for key in m_length.keys():
			mfile.write("%d\t%d\n" % (key,m_length[key]))
	with open("../data/in_length.txt","w") as ifile:
		for key in i_length.keys():
			ifile.write("%d\t%d\n" % (key,i_length[key]))
        with open("../data/out_length.txt","w") as ofile:
		for key in o_length.keys():
			ofile.write("%d\t%d\n" % (key,o_length[key]))

