#parse the given training set protein sequence file and generate two list of protein sequence and corresponding feature of i,m,o

import sys, math, argparse
from itertools import groupby


#read test fasta protein sequence file
#def read_protein(fasta):
#	fasta_list = []
#	for line in fasta:
#		if line[0]=='>':
#			try:
#				fasta_list.append(current_pro)
#			except UnboundLocalError:
#				pass
#			current_pro = [line.lstrip('>').rstrip('\n'),'']
#			print current_pro
#		else:
#			current_pro[1]+= ''.join(line.split())
#			print current_pro
#		fasta_list.append(current_pro)
#	return fasta_list


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

#parse the training protein sequence and features
def parse_training(training_file):
	seq_set=[]
	feature_set=[]
	for line in training_file:
		tab = line.strip().split('\t')
		seq_set.append(tab[0].upper())
		feature_set.append(tab[1].upper())
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
		lengthFreq[key]= float(length_set[key]/domain_total)
	return lengthFreq

# calculate the emission probability from each state to amino acids
def EmissionProb(seq_set, feature_set):
	element = len(seq_set)
	emit_dict = {}
	for x in range(element):
		size = len(seq_set[x])
		for y in range(size):
			curr_emit = feature_set[x][y]+'->'+seq_set[x][y]
			if curr_emit not in emit_dict:
				emit_dict[curr_emit] = 1
			else:
				emit_dict[curr_emit] += 1
	mem, inner, outer, non = 0, 0, 0, 0
	for key in emit_dict:
		if 'M->' in key:
			mem += emit_dict[key]
		elif 'I->' in key:
			inner += emit_dict[key]
		elif 'O->' in key:
			outer += emit_dict[key]
		else:
			non += emit_dict[key]
	for key in emit_dict:
		if 'M->' in key:
			emit_dict[key] = float(emit_dict[key])/mem
		elif 'I->' in key:
			emit_dict[key] = float(emit_dict[key])/inner
		elif 'O->' in key:
			emit_dict[key] = float(emit_dict[key])/outer
		else:
			emit_dict[key] = float(emit_dict[key])/non
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
		elif 'I->' in key:
			inner += trans_dict[key]
		else:
			outer += trans_dict[key]
	for key in trans_dict:
		if 'M->' in key:
			trans_dict[key] = float(trans_dict[key])/mem
		elif 'I->' in key:
			trans_dict[key] = float(trans_dict[key])/inner
		else:
			trans_dict[key] = float(trans_dict[key])/outer
	return trans_dict

# calculate the probability of initial state
def InitState(feature_set):
	element = len(seq_set)
	init_state = {}
	for x in range(element):
		curr_init = feature_set[x][0]
		if '.' in curr_init: continue
		if curr_init not in init_state:
			init_state[curr_init] = 1
		else:
			init_state[curr_init] += 1
	mem, inner, outer = 0, 0, 0
	for key in init_state:
		if 'M->' in key:
			mem += init_state[key]
		elif 'I->' in key:
			inner += init_state[key]
		else:
			outer += init_state[key]
	for key in init_state:
		if 'M->' in key:
			init_state[key] = float(init_state[key])/mem
		elif 'I->' in key:
			init_state[key] = float(init_state[key])/inner
		else:
			init_state[key] = float(init_state[key])/outer
	return init_state

def NullModel(seq,emit_prob):
	prob = 0.0
	for x in range(len(seq)):
		curr_emit = '.->'+seq[x]
		prob += math.log(emit_prob[curr_emit])
	return init_state

def NullModel(seq,emit_dict):
	prob = 0.0
	for x in range(len(seq)):
		curr_emit = '.->'+seq[x]
		prob += math.log(emit_dict[curr_emit])
	return prob

#calculate the maximum probability of hidden states sequence
def max_hidden(seq,len_table,emit_dict,trans_dict,initial):
	l=len(seq)
	m=len(len_table)
	print("length of seq is %d" %(l))
	states=['M','I','O']
	pro_table = [[1 for x in range(l)]for x in range(m)]#pro_table is the max_hidden states probability table
	pos_table = [[1 for x in range(l)]for x in range(m)]#pos_table is the path table for max_hidden states
	pro_table[0][0] = 0.000001 #actually 0 but since we use log, we assign this very small value for staring M
	pro_table[1][0] = initial['I']
	pro_table[2][0] = initial['O']
	for i in range(1,l):
		for j in range(m):
			'''pro_oneseq is the probability that observed segment ends in state j is continous '''
			pro_oneseg = math.log(len_table[j][i+1])+math.log(pro_table[j][0])
			for aa in seq[:i+1]:
				emit = states[j]+'->'+ aa.upper()
				pro_oneseg += math.log(emit_dict[emit])
			pro_table[j][i] = pro_oneseg
			pos_table[j][i] = (j,0)
			for k in range(1,i):
				for p in range(m):
					tmp = 0
					if p!=j:
						seg = seq[k+1:i+1]
						trans = states[p] +'->'+ states[j]
						if trans not in trans_dict.keys():
							trans_dict[trans] = 0.000001#len_table is the length distribution table
						tmp += pro_table[p][k] + math.log(trans_dict[trans]) + math.log(len_table[j][i-k])
						for aa in seg:
							emit = states[p] + '->' + aa.upper()
							tmp += math.log(emit_dict[emit])
						if pro_table[j][i]< tmp:
							pro_table[j][i] = tmp
							pos_table[j][i] = (p,k)
	print("The length of output table is %d" %(len(pos_table[0])))
	return pro_table,pos_table

def TraceBack(seq,emit_dict,pro_table,pos_table):
	x = len(seq)-1
	HiddenStates = ''
	NullScore = NullModel(seq,emit_dict)
	TMScore = max(pro_table[0][x],pro_table[1][x],pro_table[2][x])
	if NullScore > TMScore:
		HiddenStates = '.'*len(seq)
		return HiddenStates

	if max(pro_table[0][x],pro_table[1][x],pro_table[2][x]) == pro_table[0][x]:
		record = pos_table[0][x]
		state = 'M'
	elif max(pro_table[0][x],pro_table[1][x],pro_table[2][x]) == pro_table[1][x]:
		record = pos_table[1][x]
		state = 'I'
	else:
		record = pos_table[2][x]
		state = 'O'
	repeats = x - record[1]
	HiddenStates += state*repeats
	pre_state = record[0]
	pre_pos = record[1]

	while True:
		record = pos_table[pre_state][pre_pos]
		if pre_state == 0:
			state = 'M'
		elif pre_state == 1:
			state = 'I'
		else:
			state = 'O'
		repeats = pre_pos - record[1]
		HiddenStates += state*repeats
		pre_state = record[0]
		pre_pos = record[1]
		if pre_pos == 0:
			HiddenStates += state
			break

	HiddenStates = HiddenStates[::-1]
	return HiddenStates


###======read traing file and generate GMMM model ======
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--fasta_file',required=True)
	parser.add_argument('-o','--out_file',required=True)
	parser.add_argument('-l','--length_table')
	args = parser.parse_args()

	with open("../data/TMseq.ffa","r") as training_file:
		seq_set,feature_set = parse_training(training_file)
		m_length,i_length,o_length = generate_length(feature_set)
		emit_prob = EmissionProb(seq_set, feature_set)
		trans_prob = TransitionProb(feature_set)
		init_state = InitState(feature_set)
#	print(freq_table)
	#with open(args.fasta_file,'r') as test:
		#### test_seq_labels contains the line labels
	test_class = classFASTA(args.fasta_file)
	test_seq = [ x[1] for x in test_class.readFASTA() ]
	test_labels  = [ x[0] for x in test_class.readFASTA() ]


#add pesudocount to length distribution for each domain
	for i in range(len(test_seq)):
		if i+1 not in m_length.keys():
			m_length[i+1] = 0.01
		else:
			m_length[i+1] += 0.01
		if i+1 not in i_length.keys():
			i_length[i+1] = 0.01
		else:
			i_length[i+1] +=0.01
		if i+1 not in o_length.keys():
			o_length[i+1] = 0.01
		else:
			o_length[i+1] += 0.01

	m_len = lengthFrequency(m_length)
	i_len = lengthFrequency(i_length)
	o_len = lengthFrequency(o_length)
	len_table = [m_len,i_len,o_len]
	maxpro_table,pos_table = max_hidden(test_seq,len_table,emit_prob,trans_prob,init_state)

	hidden_states = TraceBack(test_seq,emit_prob,maxpro_table,pos_table)

	with open("../data/maxpro_table","w")as ofile:
		for i in range(len(maxpro_table)):
			for j in range(len(maxpro_table[i])):
				ofile.write(str(maxpro_table[i][j])+'\t')
			ofile.write('\n')
	with open("../data/pos_table","w")as ofile:
		for i in range(len(pos_table)):
			for j in range(len(pos_table[i])):
				ofile.write(str(pos_table[i][j])+'\t')
			ofile.write('\n')


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
		for key in sorted(emit_prob.iterkeys()):
			ofile.write("%s\t%f\n" % (key,emit_prob[key]))
	with open("../data/trans_probabity.txt","w") as ofile:
		for key in sorted(trans_prob.iterkeys()):
			ofile.write("%s\t%f\n" % (key,trans_prob[key]))
	with open("../data/initial_state.txt","w") as ofile:
		for key in sorted(init_state.iterkeys()):
			ofile.write("%s\t%f\n" % (key,init_state[key]))
	with open(args.out_file,"w") as ofile:
		ofile.write("%s\n" % (test_seq))
		ofile.write("%s\n" % (hidden_states))
