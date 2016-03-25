#this script is the implementation of pruning algorithm
from __future__ import division
import re
import os

#default setting and data for this project
base_length = 0.1
base_matrix=[[0.9,0.05,0.025,0.025],[0.05,0.9,0.025,0.025],[0.025,0.025,0.9,0.05],[0.025,0.025,0.05,0.9]]#the substitution matrix for length = 0.1
background_freq = [0.25,0.25,0.25,0.25]#assume equal probabilty for A,C,G,T

# This function reads in a Newick formatted file

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'


def read_newick(newick_list):
	'''
	this function takes a list of newick lines as an arguement
	and returns each item in the list as a nested tuple
	The function assumes you have a sorted Newick object
	'''
	for nwk in newick_list:
		nested_levels = nwk.count(')') + 1
		nwk_list =  re.split(r'[,:]+', nwk)
		#print nwk_list
		nwk_list = [line.translate(None, '()') for line in nwk_list]
		number_comparisons = sum(range(0, nested_levels, 1))
		return_list = []
		for x in range(0,number_comparisons):
			if x == 0:
				label = str(nwk_list[x]) + '_' + str(nwk_list[x+2])
				return_tuple = (label, float(nwk_list[1]), float(nwk_list[3]))
				return_list.append(return_tuple)
			elif x == 1 :
				label = str(nwk_list[x+1]) + '_' + str(nwk_list[x+4])
				return_tuple = (label, float(nwk_list[1]) + float(nwk_list[4]) , float(nwk_list[6]))
				return_list.append(return_tuple)
			elif x == 2 :
				label = str(nwk_list[x]) + '_' + str(nwk_list[x+3])
				return_tuple = (label, float(nwk_list[3]) + float(nwk_list[4]) , float(nwk_list[6]))
				return_list.append(return_tuple)

		print return_list
		return return_list


#this function does 4x4 matrix multiplication
def matrix_multiply(matrixA,matrixB):
	output = [[0 for i in range(4)]for i in range(4)]
	for i in range(4):
		for j in range(4):
			output[i][j] = matrixA[i][0]*matrixB[0][j]+matrixA[i][1]*matrixB[1][j]+matrixA[i][2]*matrixB[2][j]+matrixA[i][3]*matrixB[3][j]

	return output

#this function take input of base_substitution matrix and the value of t, outputp(t)
def Compute_substitution(base_matrix,length):
	relative_length = int(length /0.1)
	'''initial of the output matrix in length t '''
	new_matrix =[[0 for i in range(4)]for i in range(4)]
	if relative_length == 1:
		return base_matrix
	elif relative_length == 2:
		return matrix_multiply(base_matrix,base_matrix)
	else:
		new_matrix = matrix_multiply(base_matrix,base_matrix)
		for i in range(relative_length-2):
			new_matrix = matrix_multiply(new_matrix,base_matrix)
		return new_matrix

#function of Pruning Algorithm to calculate the A,C,G,T freq for the parental node from its two children
def PruningAlgrtm(node1,node2,matrix1,matrix2):
	# node1, node2 are two children of one node
	print node1, node2
	node1_A, node1_C, node1_G, node1_T = node1[0], node1[1], node1[2], node1[3]
	node2_A, node2_C, node2_G, node2_T = node2[0], node2[1], node2[2], node2[3]
	#print node2_A
	# freq of parental node is A
	nodeP_A = (node1_A*matrix1[0][0]+node1_C*matrix1[0][1]+node1_G*matrix1[0][2]+node1_T*matrix1[0][3]) * (node2_A*matrix2[0][0]+node2_C*matrix2[0][1]+node2_G*matrix2[0][2]+node2_T*matrix2[0][3])
	# freq of parental node is C
	nodeP_C = (node1_A*matrix1[1][0]+node1_C*matrix1[1][1]+node1_G*matrix1[1][2]+node1_T*matrix1[1][3]) * (node2_A*matrix2[1][0]+node2_C*matrix2[1][1]+node2_G*matrix2[1][2]+node2_T*matrix2[1][3])
	# freq of parental node is G
	nodeP_G = (node1_A*matrix1[2][0]+node1_C*matrix1[2][1]+node1_G*matrix1[2][2]+node1_T*matrix1[2][3]) * (node2_A*matrix2[2][0]+node2_C*matrix2[2][1]+node2_G*matrix2[2][2]+node2_T*matrix2[2][3])
	# freq of parental node is T
	nodeP_T = (node1_A*matrix1[3][0]+node1_C*matrix1[3][1]+node1_G*matrix1[3][2]+node1_T*matrix1[3][3]) * (node2_A*matrix2[3][0]+node2_C*matrix2[3][1]+node2_G*matrix2[3][2]+node2_T*matrix2[3][3])
	nodeP = [nodeP_A,nodeP_C,nodeP_G,nodeP_T]
	#print nodeP
	return nodeP

#A function to test the accuray of the program using class example
def TestExmpl():
	#test_tree = '(((A:0.2, B:0.2):0.1, C:0.2):0.1, (D:0.2, E:0.2):0.1)'
	#test_leaves = ['T', 'C', 'A', 'C', 'C']
	# matrix order is based on 'A, C, G, T'   ('T, C, A, G' shown in the class slide)
	pass
test_matrix_01= [[0.906563,0.023791,0.045855,0.023791],[0.023791,0.906563,0.023791,0.045855],[0.045855,0.023791,0.906563,0.023791],[0.023791,0.045855,0.023791,0.906563]]
test_matrix_02 = Compute_substitution(test_matrix_01,0.2)

#this function read the tree_file in Newick Standard format, then output a binary tree in nested lists
def read_tree(tree_file):
	with open(tree_file,'r') as ifile:
		tree_line = ifile.read().strip()



#test of Compute_substitution function

	leaf1 = [0,0,0,1]
	leaf2 = [0,1,0,0]
	leaf3 = [1,0,0,0]
	leaf4 = [0,1,0,0]
	leaf5 = [0,1,0,0]


	# length of branch to node to 6,7,8 is 0.1, others are 0.2
	node7 = PruningAlgrtm(leaf1,leaf2,test_matrix_02,test_matrix_02)
	node6 = PruningAlgrtm(node7,leaf3,test_matrix_01,test_matrix_02)
	node8 = PruningAlgrtm(leaf4,leaf5,test_matrix_02,test_matrix_02)
	root = PruningAlgrtm(node6,node8,test_matrix_01,test_matrix_01)
	return root




###### use the program to work on project data #####
if __name__ == "__main__":
	#print TestExmpl()  #test the accuray of the program

	'''
	pattern = re.compile(r"\b[0-9]+(?:\.[0-9]+)?\b")
	branch_lengths = pattern.findall(tree)

	#print(sum(test_matrix_01[0]),sum(test_matrix_01[1]),sum(test_matrix_01[2]),sum(test_matrix_01[3])) #test the matrix
	'''


#new_matrix = Compute_substitution(test_matrix,0.3)






#print(new_matrix)
#print(sum(new_matrix[0]))
print read_tree("../tree_file")

# Testing newick reader

#with open(mydir + 'tree_file.txt') as IN:
#	read_tree(IN)
	#trees = [line.split('\n') for line in IN.read().strip().split('\n\n')]
#trees_test = read_newick(trees[0])
