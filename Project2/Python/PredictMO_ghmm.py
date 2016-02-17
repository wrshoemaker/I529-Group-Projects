#parse the given training set protein sequence file and generate two list of protein sequence and corresponding feature of i,m,o

import sys
from itertools import groupby
import math
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

# calculate the emission probability from each state to amino acids
def EmissionProb(seq_set, feature_set):
	element = len(seq_set)
	emit_dict = {}
	for x in range(element):	
		size = len(seq_set[x])
		for y in range(size):
			curr_emit = feature_set[x][y]+'->'+seq_set[x][y]
			if '.' in curr_emit: continue
			if curr_emit not in emit_dict:
				emit_dict[curr_emit] = 1
			else:
				emit_dict[curr_emit] += 1
	mem, inner, outer = 0, 0, 0
	for key in emit_dict:
		if 'M->' in key:
			mem += emit_dict[key]
		elif 'i->' in key:
			inner += emit_dict[key]
		else:
			outer += emit_dict[key]
	for key in emit_dict:
		if 'M->' in key:
			emit_dict[key] = float(emit_dict[key])/mem
		elif 'i->' in key:
			emit_dict[key] = float(emit_dict[key])/inner
		else:
			emit_dict[key] = float(emit_dict[key])/outer
	return emit_dict

# calculate the transition probability between states
def TransitionProb(feature_set):
	element = len(seq_set)
	trans_dict = {}
	for x in range(element):	
		size = len(seq_set[x])
		for y in range(size-1):
			if feature_set[x][y] == feature_set[x][y+1]: continue
			curr_trans = feature_set[x][y]+'->'+feature_set[x][y+1]
			if '.' in curr_trans: continue
			if curr_trans not in trans_dict:
				trans_dict[curr_trans] = 1
			else:
				trans_dict[curr_trans] += 1
	mem, inner, outer = 0, 0, 0
	for key in trans_dict:
		if 'M->' in key:
			mem += trans_dict[key]
		elif 'i->' in key:
			inner += trans_dict[key]
		else:
			outer += trans_dict[key]
	for key in trans_dict:
		if 'M->' in key:
			trans_dict[key] = float(trans_dict[key])/mem
		elif 'i->' in key:
			trans_dict[key] = float(trans_dict[key])/inner
		else:
			trans_dict[key] = float(trans_dict[key])/outer
	return trans_dict
#calcualte the max
def max_hidden(seq,trans_dict,emit_dict,freq_table):
	l=len(seq)

		
###======read traing file and generate GMMM model ======
if __name__ == "__main__":
	with open("../data/TMseq.ffa","r") as training_file:
		seq_set,feature_set = parse_training(training_file)
		m_length,i_length,o_length = generate_length(feature_set)
		m_freq = lengthFrequency(m_length)
		i_freq = lengthFrequency(i_length)
		o_freq = lengthFrequency(o_length)
		emit_prob = EmissionProb(seq_set, feature_set)
		trans_prob = TransitionProb(feature_set)

	with open("../data/mem_length.txt","w") as mfile:
		for key in m_length.keys():
			mfile.write("%d\t%d\n" % (key,m_length[key]))
	with open("../data/in_length.txt","w") as ifile:
		for key in i_length.keys():
			ifile.write("%d\t%d\n" % (key,i_length[key]))
	with open("../data/out_length.txt","w") as ofile:
		for key in o_length.keys():
			ofile.write("%d\t%d\n" % (key,o_length[key]))
	with open("../data/emit_probabity.txt","w") as ofile:
		for key in emit_prob.keys():
			ofile.write("%s\t%f\n" % (key,emit_prob[key]))
	with open("../data/trans_probabity.txt","w") as ofile:
		for key in trans_prob.keys():
			ofile.write("%s\t%f\n" % (key,trans_prob[key]))

