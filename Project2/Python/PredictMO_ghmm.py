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

#this class check input format and parse fasta file into two list containing index and sequence
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

#parse the training protein into peptide sequence and features
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
	m_length={}#key is length, value is frequency
	i_length={}
	o_length={}
	for sequence in feature_set:
		tmp=[''.join(g) for k,g in groupby(sequence)] #parse the feature sequence into tmp=['mmm','iii','oooo'..]
		for each in tmp:
			if 'm' in each.lower():#the membrane segment
				if len(each) not in m_length.keys():
					m_length[len(each)] = 1
				else:
					m_length[len(each)] += 1
			elif 'i' in each.lower():#inside segment
				if len(each) not in i_length.keys():
					i_length[len(each)] = 1
				else:
					i_length[len(each)] += 1
			elif 'o' in each.lower():#outside segment
				if len(each) not in o_length.keys():
					o_length[len(each)] = 1
				else:
					o_length[len(each)] += 1
	#return three dictionaries
	return m_length,i_length,o_length

#calculate the frequency of length for i,m,o domain
def lengthFrequency(length_set):
	'''the length_set is the length distribution dictionary for each part, key is length, value is frequency
	like m_length[2] = 10 means a segment 'mm' occurs 10 times'''
	lengthFreq={}
	domain_total=sum(length_set.values())
	for key in length_set:
		lengthFreq[key]= float(length_set[key])/domain_total
	return lengthFreq

# calculate the emission probability from each state to amino acids
def EmissionProb(seq_set, feature_set):
	'''seq_set is a list of peptide sequences, feature_set is a list of features'''
	element = len(seq_set)#number of protein sequences
	emit_dict = {}#the emission table in format emit_dict['M->A']=
	for x in range(element):
		size = len(seq_set[x])
		for y in range(size):
			curr_emit = feature_set[x][y]+'->'+seq_set[x][y]
			if curr_emit not in emit_dict:
				emit_dict[curr_emit] = 1
			else:
				emit_dict[curr_emit] += 1
	mem, inner, outer, non = 0, 0, 0, 0#sum for each part
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
def TransitionProb():#four possible states O->MO->I->MI->O and only four possible transition
	trans_dict={}
	states=['O','MO','I','MI'] #list of four states
	for start in states:
		for end in states:
			curr_trans = start + '->' + end
			trans_dict[curr_trans] = 0.00000000001#give a very small pseudocount for all transition
	'''actually only four transition with probabilty 1 are allowed'''
	trans_dict['O->MO'] = 1
	trans_dict['MO->I'] = 1
	trans_dict['I->MI'] = 1
	trans_dict['MI->O'] = 1
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
	l=len(seq)#length for each test sequence
	m=len(len_table)#number of hidden states,len_table is list of three dictionaries 
	states=['MI','MO','I','O']
	pro_table = [[1 for x in range(l)]for x in range(m)]#pro_table is the max_hidden states probability table
	pos_table = [[1 for x in range(l)]for x in range(m)]#pos_table is the path table for max_hidden states
	pro_table[0][0] = 0.00000000001 #actually 0 but since we use log, we assign this very small value for starting M
	pro_table[1][0] = 0.00000000001
	pro_table[2][0] = initial['I']
	pro_table[3][0] = initial['O']
	for i in range(1,l):
		for j in range(m):
			'''pro_oneseq is the probability that observed segment ends in state j is continous '''
			pro_oneseg = math.log(len_table[j][i+1])+math.log(pro_table[j][0])
			for aa in seq[:i+1]:
				if states[j]=='MI' or states[j]=='MO':#'MI'and'MO' have same emission probabilty which is emit_dict['M->?']
					emit = 'M->'+ aa.upper()
				else:
					emit = states[j]+'->'+ aa.upper()
				pro_oneseg += math.log(emit_dict[emit])
			pro_table[j][i] = float(pro_oneseg)#assign the probabilty of first conditon to [j][i] at first
			pos_table[j][i] = (j,0)#ths position for first conditon is to start from position 0
			for k in range(1,i):
				for p in range(m):
					tmp = 0#the probability for second condition
					if p!=j:#the second condition is to start from any cells not in the same row with final state
						seg = seq[k+1:i+1]#the duration sequence in state j
						trans = states[p] +'->'+ states[j]#state transition from p to j
						tmp += pro_table[p][k] + math.log(trans_dict[trans]) + math.log(len_table[j][i-k])#the preceding state [p][k] and transition probability and duration probability
						for aa in seg:#the emission of aa in state j, I ignored a typo here that I wrote k instead of j
							if j==0 or j==1:#combine emission of MI,MO into M
								emit = 'M->'+ aa.upper()
							else:
								emit = states[j] + '->' + aa.upper()
							tmp += math.log(emit_dict[emit])
						if pro_table[j][i]< tmp:
							pro_table[j][i] = tmp #tmp max value
							pos_table[j][i] = (p,k)#postion of preceding state
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

	if max(pro_table[0][x],pro_table[1][x],pro_table[2][x],pro_table[3][x]) == pro_table[0][x]:
		record = pos_table[0][x]
		state = 'M'
	elif max(pro_table[0][x],pro_table[1][x],pro_table[2][x],pro_table[3][x]) == pro_table[1][x]:
		record = pos_table[1][x]
		state = 'M'
	elif max(pro_table[0][x],pro_table[1][x],pro_table[2][x],pro_table[3][x]) == pro_table[2][x]:
		record = pos_table[2][x]
		state = 'I'
	else:
		record = pos_table[3][x]
		state = 'O'
	repeats = x - record[1]
	HiddenStates += state*repeats
	pre_state = record[0]
	pre_pos = record[1]

	while True:
		record = pos_table[pre_state][pre_pos]
		if pre_state == 0 or pre_state == 1:
			state = 'M'
		elif pre_state == 2:
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
		m_freq = lengthFrequency(m_length)
		i_freq = lengthFrequency(i_length)
		o_freq = lengthFrequency(o_length)
#	with open("../data/mem_freq.txt","w") as ofile:
#		for key in m_freq.keys():
#			ofile.write(str(key)+"\t"+str(m_freq[key])+"\n")
#	with open("../data/in_freq.txt","w") as ofile:
#		for key in i_freq.keys():
#			ofile.write(str(key)+"\t"+str(i_freq[key])+"\n")
#	with open("../data/out_freq.txt","w") as ofile:
#		for key in o_freq.keys():
#			ofile.write(str(key)+"\t"+str(o_freq[key])+"\n")

		emit_prob = EmissionProb(seq_set, feature_set)
		trans_prob = TransitionProb()
		init_state = InitState(feature_set)
#	print(freq_table)
	#with open(args.fasta_file,'r') as test:
		#### test_seq_labels contains the line labels
	test_class = classFASTA(args.fasta_file)
	test_seqs = [ x[1] for x in test_class.readFASTA() ]
	test_labels  = [ x[0] for x in test_class.readFASTA() ]

	index = 0

	for test_seq in test_seqs:#test_seq is each test protein sequence
		tmpm_length,tmpi_length,tmpo_length = {}, {}, {}#three dictionaries for m,i,o length distribution
		for key in m_length.keys():#use training prior knowledge
			tmpm_length[key] = m_length[key]
		for key in i_length.keys():
			tmpi_length[key] = i_length[key]
		for key in o_length.keys():
			tmpo_length[key] = o_length[key]
#		print("the length of original mem")
#		print(tmpm_length)
		test_id = test_labels[index]
		index += 1
		print("the length of this testing sequence is %d" %(len(test_seq)))
		#add pesudocount to length distribution for each domain
		for i in range(len(test_seq)):
			if i+1 not in tmpm_length.keys():
				tmpm_length[i+1] = 0.00000000001
			else:
				tmpm_length[i+1] += 0.00000000001
			if i+1 not in tmpi_length.keys():
				tmpi_length[i+1] = 0.00000000001
			else:
				tmpi_length[i+1] += 0.00000000001
			if i+1 not in tmpo_length.keys():
				tmpo_length[i+1] = 0.00000000001
			else:
				tmpo_length[i+1] += 0.00000000001

		tmpm_len = lengthFrequency(tmpm_length)#calculate probability table for each length table
		tmpi_len = lengthFrequency(tmpi_length)
		tmpo_len = lengthFrequency(tmpo_length)
#check length probabilty
		for check in [tmpm_len,tmpi_len,tmpo_len]:#check if sum of probability is 1
			print("the sum of probabilty is %f" %(sum(check.values())))
		len_table = [tmpm_len,tmpm_len,tmpi_len,tmpo_len]#list of length dictionaries in order MI,MO,I,O. MI and Mo are actually same
		#for domain in tmpm_len.keys():
		#	print("probability of fragment in length %d in mem is %f" %(domain,tmpm_len[domain]))
		#for domain in tmpi_len.keys():
			#print("probability of fragment in length %d in in is %f"%(domain,tmpi_len[domain]))
		#for domain in tmpo_len.keys():
			#print("probability of fragment in length %d in out is %f"%(domain,tmpo_len[domain]))
#calculate the viterbi value table and path table,then produce string of hidden_states
		maxpro_table,pos_table = max_hidden(test_seq,len_table,emit_prob,trans_prob,init_state)
		hidden_states = TraceBack(test_seq,emit_prob,maxpro_table,pos_table)
#append the hidden_states to output file following seq index and protein sequence
		with open(args.out_file,"a") as ofile:
			ofile.write(">%s\n" % (test_id))
			ofile.write("%s\n" % (test_seq))
			ofile.write("%s\n" % (hidden_states))