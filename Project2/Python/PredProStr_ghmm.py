#parse the given training set protein sequence file and generate two list of protein sequence and corresponding feature of i,m,o

import sys, math, argparse
from itertools import groupby
import re
		
#this class check input format and parse fasta file into two list containing index and sequence
class classFASTA:

	def __init__(self, fileFASTA):
		self.fileFASTA = fileFASTA

	def readFASTA(self):
		return self.ParseFASTA(self.fileFASTA)

	def ParseFASTA(self, fileFASTA):
		fasta_list=[]
		for line in fileFASTA:
			line = line.rstrip()
			if '<>' in line:
				curr_data = [[],[]]
			elif 'end' in line:
				curr_data[0] = ''.join(curr_data[0])
				curr_data[1] = ''.join(curr_data[1])
				fasta_list.append(curr_data)
			elif (re.search('[a-zA-Z]', line)):	
				curr_data[0].append(line.split()[0])
				curr_data[1].append(line.split()[1])
		return fasta_list

#table of peptide length frequency in both membrane and non-membrane domain
def generate_length(feature_set):
	h_length={}#key is length, value is frequency
	e_length={}
	c_length={}
	for sequence in feature_set:
		tmp=[''.join(g) for k,g in groupby(sequence)] #parse the feature sequence into tmp=['mmm','iii','oooo'..]
		for each in tmp:
			if 'h' in each.lower():#the membrane segment
				if len(each) not in h_length.keys():
					h_length[len(each)] = 1
				else:
					h_length[len(each)] += 1
			elif 'e' in each.lower():#inside segment
				if len(each) not in e_length.keys():
					e_length[len(each)] = 1
				else:
					e_length[len(each)] += 1
			elif '_' in each.lower():#outside segment
				if len(each) not in c_length.keys():
					c_length[len(each)] = 1
				else:
					c_length[len(each)] += 1
	#return three dictionaries
	return h_length,e_length,c_length

#calculate the frequency of length for i,m,o domain
def lengthFrequency(length_set):
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
	alpha, beta, coli = 0, 0, 0
	for key in emit_dict:
		if 'h->' in key:
			alpha += emit_dict[key]
		elif 'e->' in key:
			beta += emit_dict[key]
		elif '_->' in key:
			coli += emit_dict[key]
	for key in emit_dict:
		if 'h->' in key:
			emit_dict[key] = float(emit_dict[key])/alpha
		elif 'e->' in key:
			emit_dict[key] = float(emit_dict[key])/beta
		elif '_->' in key:
			emit_dict[key] = float(emit_dict[key])/coli
	return emit_dict

# calculate the transition probability between states
def TransitionProb(feature):
	trans_dict={}
	states=['h','e','_'] 
	trans_dict = {'h->e':0, 'h->_':0, 'e->h':0, 'e->_':0, '_->h':0, '_->e':0}
	for x in range(len(feature)):
		seq = feature[x]
		for y in range(len(seq)-1):
			if seq[y] != seq[y+1]:
				curr_trans = seq[y]+'->'+seq[y+1]
				for key in trans_dict:
					if key == curr_trans:
						trans_dict[curr_trans] += 1
	alpha = trans_dict['h->e']+trans_dict['h->_']
	beta = trans_dict['e->h']+trans_dict['e->_']
	coli = trans_dict['_->h']+trans_dict['_->e']
	trans_dict['h->e'] = float(trans_dict['h->e'])/alpha
	trans_dict['h->_'] = float(trans_dict['h->_'])/alpha
	trans_dict['e->h'] = float(trans_dict['e->h'])/beta
	trans_dict['e->_'] = float(trans_dict['e->_'])/beta
	trans_dict['_->h'] = float(trans_dict['_->h'])/coli
	trans_dict['_->e'] = float(trans_dict['_->e'])/coli
	return trans_dict

# calculate the probability of initial state
def InitState(feature_set):
	element = len(seq_set)
	init_state = {}
	for x in range(element):
		curr_init = feature_set[x][0]
		if curr_init not in init_state:
			init_state[curr_init] = 1
		else:
			init_state[curr_init] += 1
	alpha, beta, coli = 0, 0, 0
	for key in init_state:
		if 'h' in key:
			alpha += init_state[key]
		elif 'e' in key:
			beta += init_state[key]
		else:
			coli += init_state[key]
	for key in init_state:
		if 'h' in key:
			init_state[key] = float(init_state[key])/alpha
		elif 'e' in key:
			init_state[key] = float(init_state[key])/beta
		else:
			init_state[key] = float(init_state[key])/coli
	return init_state


#calculate the maximum probability of hidden states sequence
def max_hidden(seq,len_table,emit_dict,trans_dict,initial):
	l=len(seq)#length for each test sequence
	m=len(len_table)#number of hidden states,len_table is list of three dictionaries 
	states=['h','e','_']
	pro_table = [[1 for x in range(l)]for x in range(m)]#pro_table is the max_hidden states probability table
	pos_table = [[1 for x in range(l)]for x in range(m)]#pos_table is the path table for max_hidden states
	pro_table[0][0] = 0.00000000001 #actually 0 but since we use log, we assign this very small value for starting M
	pro_table[1][0] = 0.00000000001
	pro_table[2][0] = initial['_']
	for i in range(1,l):
		for j in range(m):
			'''pro_oneseq is the probability that observed segment ends in state j is continous '''
			pro_oneseg = math.log(len_table[j][i+1])+math.log(pro_table[j][0])  ##initial probabilty * length ditribution probability (from position 0)
			for aa in seq[:i+1]:
				emit = states[j]+'->'+ aa.upper()
				pro_oneseg += math.log(emit_dict[emit])
			pro_table[j][i] = float(pro_oneseg)  ##assign the probabilty of first conditon to [j][i] at first
			pos_table[j][i] = (j,0)  ##ths position for first conditon is to start from position 0
			for k in range(1,i):
				for p in range(m):
					tmp = 0  ##the probability for second condition
					if p!=j:  ##the second condition is to start from any cells not in the same row with final state
						seg = seq[k+1:i+1]#the duration sequence in state j
						trans = states[p] +'->'+ states[j]#state transition from p to j
						tmp += pro_table[p][k] + math.log(trans_dict[trans]) + math.log(len_table[j][i-k])#the preceding state [p][k] and transition probability and duration probability
						for aa in seg:#the emission of aa in state j, I ignored a typo here that I wrote k instead of j
							emit = states[j] + '->' + aa.upper()
							tmp += math.log(emit_dict[emit])
						if pro_table[j][i]< tmp:
							pro_table[j][i] = tmp #tmp max value
							pos_table[j][i] = (p,k)#postion of preceding state
	return pro_table,pos_table

def TraceBack(seq,emit_dict,pro_table,pos_table):
	x = len(seq)-1
	HiddenStates = ''
	TMScore = max(pro_table[0][x],pro_table[1][x],pro_table[2][x])

	if max(pro_table[0][x],pro_table[1][x],pro_table[2][x]) == pro_table[0][x]:
		record = pos_table[0][x]
		state = 'h'
	elif max(pro_table[0][x],pro_table[1][x],pro_table[2][x]) == pro_table[1][x]:
		record = pos_table[1][x]
		state = 'e'
	elif max(pro_table[0][x],pro_table[1][x],pro_table[2][x]) == pro_table[2][x]:
		record = pos_table[2][x]
		state = '_'
	repeats = x - record[1]
	HiddenStates += state*repeats
	pre_state = record[0]
	pre_pos = record[1]

	while True:
		record = pos_table[pre_state][pre_pos]
		if pre_state == 0:
			state = 'h'
		elif pre_state == 1:
			state = 'e'
		elif pre_state == 2:
			state = '_'
		repeats = pre_pos - record[1]
		HiddenStates += state*repeats
		pre_state = record[0]
		pre_pos = record[1]
		if pre_pos == 0:
			HiddenStates += state
			break

	HiddenStates = HiddenStates[::-1]
	return HiddenStates

def AccuracyTest(real_states, pred_states):
	total, positive = 0, 0
	for x in range(len(real_states)):
		total += 1
		if pred_states[x] == real_states[x]:
			positive += 1
	accuracy = str(float(positive)/total*100)+'%'
	return accuracy

###======read traing file and generate GMMM model ======
if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--fasta_file',required=True)
	parser.add_argument('-o','--out_file',required=True)
	parser.add_argument('-l','--length_table')
	args = parser.parse_args()
	
	with open("training.txt","r") as training_file:
		train_set = classFASTA(training_file)
		seq_set = [ x[0] for x in train_set.readFASTA() ]
		
	with open("training.txt","r") as training_file:
		train_set = classFASTA(training_file)
		feature_set = [ x[1] for x in train_set.readFASTA() ]

	h_length,e_length,c_length = generate_length(feature_set)
	h_freq = lengthFrequency(h_length)
	e_freq = lengthFrequency(e_length)
	c_freq = lengthFrequency(c_length)

	with open("alpha_freq.txt","w") as ofile:
		for key in h_freq.keys():
			ofile.write(str(key)+"\t"+str(h_freq[key])+"\n")
	with open("beta_freq.txt","w") as ofile:
		for key in e_freq.keys():
			ofile.write(str(key)+"\t"+str(e_freq[key])+"\n")
	with open("coli_freq.txt","w") as ofile:
		for key in c_freq.keys():
			ofile.write(str(key)+"\t"+str(c_freq[key])+"\n")

	emit_prob = EmissionProb(seq_set, feature_set)
	trans_prob = TransitionProb(feature_set)
	init_state = InitState(feature_set)
	
	test_class = classFASTA(args.fasta_file)
	
	with open(args.fasta_file,"r") as test_set:
		test_set = classFASTA(test_set)
		test_seqs = [ x[0] for x in test_set.readFASTA() ]
		
	with open(args.fasta_file,"r") as test_set:
		test_set = classFASTA(test_set)
		test_states = [ x[1] for x in test_set.readFASTA() ]	

	index = 0
	
	for i in range(len(test_seqs)):#test_seq is each test protein sequence
		test_seq = test_seqs[i]
		test_label = test_states[i]
		tmph_length,tmpe_length,tmpc_length = {}, {}, {}#three dictionaries for m,i,o length distribution
		for key in h_length.keys():#use training prior knowledge
			tmph_length[key] = h_length[key]
		for key in e_length.keys():
			tmpe_length[key] = e_length[key]
		for key in c_length.keys():
			tmpc_length[key] = c_length[key]

		for i in range(len(test_seq)):
			if i+1 not in tmph_length.keys():
				tmph_length[i+1] = 0.00000000001
			else:
				tmph_length[i+1] += 0.00000000001
			if i+1 not in tmpe_length.keys():
				tmpe_length[i+1] = 0.00000000001
			else:
				tmpe_length[i+1] += 0.00000000001
			if i+1 not in tmpc_length.keys():
				tmpc_length[i+1] = 0.00000000001
			else:
				tmpc_length[i+1] += 0.00000000001

		tmph_len = lengthFrequency(tmph_length)#calculate probability table for each length table
		tmpe_len = lengthFrequency(tmpe_length)
		tmpc_len = lengthFrequency(tmpc_length)
		len_table = [tmph_len,tmpe_len,tmpc_len]#list of length dictionaries 
		
		
		#check length probabilty
		#for check in [tmph_len,tmpe_len,tmpc_len]:#check if sum of probability is 1
		#	print("the sum of probabilty is %f" %(sum(check.values())))
		#for domain in tmpm_len.keys():
		#	print("probability of fragment in length %d in mem is %f" %(domain,tmpm_len[domain]))
		#for domain in tmpi_len.keys():
			#print("probability of fragment in length %d in in is %f"%(domain,tmpi_len[domain]))
		#for domain in tmpo_len.keys():
			#print("probability of fragment in length %d in out is %f"%(domain,tmpo_len[domain]))
		#calculate the viterbi value table and path table,then produce string of hidden_states
		
		#print len_table
		#print len(emit_prob)
		#print trans_prob
		#print init_state
		
		maxpro_table,pos_table = max_hidden(test_seq,len_table,emit_prob,trans_prob,init_state)
		pred_states = TraceBack(test_seq,emit_prob,maxpro_table,pos_table)
		
		index += 1
		real_states = ''.join(test_label)
		accuracy = AccuracyTest(real_states, pred_states)
		print ('%d\t%s') % (index, accuracy)
		
		#append the hidden_states to output file following seq index and protein sequence
		with open(args.out_file,"a") as ofile:
			ofile.write("%s\n" % (test_seq))
			ofile.write("%s\n" % (pred_states))
	