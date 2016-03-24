#this test.py is to parse Newick tree file and build a binary tree data structure 

from itertools import groupby


class Node:
	def __init__(self,length):	
		self.length = length#length to its parent
		self.left = None
		self.right = None
		self.value = None

class BinaryTree:
	def __init__(self):
		self.root = None
		self.left = None
		self.right = None
		 
class Newick:
	def __init__(self,newick):
		self.newick = newick
	def divide(self):
		if ',' not in self.newick:
			length = float(self.newick.split(':')[1])
			tmp = Node(length)
			return tmp
		else:
			tab = [''.join(g) for k,g in groupby(self.newick)]
			number_left = len(tab[0]) - 1
			number_right = 0
			total = len(self.newick)

			for i in range(total):
				if  self.newick[i] == ')':
					number_right += 1
				if number_right == number_left:
					for j in range(i+1,total):
						if self.newick[j] == ',':
							break
					break

			k = self.newick.rfind(':')
			length = float(self.newick[k+1:])
			ancestor = Node(length)
			tmp_left = self.newick[1:j]
			tmp_right = self.newick[j+1:k-1]
			left_sub = Newick(tmp_left)#left sub tree
			right_sub = Newick(tmp_right)#right sub tree
		
			ancestor.left = left_sub.divide()
			ancestor.right = right_sub.divide()
			return ancestor



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

#implementation of pruning algorithm
def pruning(base_matrix,root):
	'''base_matrix is the substitution matrix for length of 0.1 or unit length in general.
	root is the phylogeny tree and column_text is one column of the alignment '''
	if root.left == None and root.right == None:
		pruning(base_matirx,root.left)
		pruning(base_matrix,root.right)
		root.value = [0,0,0,0]
	left_matrix = Compute_substitution(base_matrix,root.left.length)
	right_matrix = Compute_substitution(base_matrix,root.right.length)
	root.value[0] = root.left.value[0]*left_matrix[0][0]+root.left.value[1]*left_matrix[1][0]+root.left.value[2]*left_matrix[2][0]+root.left.value[3]*left_matrix[3][0]
	root.value[1] = root.left.value[0]*left_matrix[0][1]+root.left.value[1]*left_matrix[1][1]+root.left.value[2]*left_matrix[2][1]+root.left.value[3]*left_matrix[3][1]
        root.value[2] = root.left.value[0]*left_matrix[0][2]+root.left.value[1]*left_matrix[1][2]+root.left.value[2]*left_matrix[2][2]+root.left.value[3]*left_matrix[3][2]
        root.value[3] = root.left.value[0]*left_matrix[0][3]+root.left.value[1]*left_matrix[1][3]+root.left.value[2]*left_matrix[2][3]+root.left.value[3]*left_matrix[3][3]
        
                
####### use the simple test Newick string to check result		
test = '((Human:0.3,Chimpanzee:0.2):0.1,Gorilla:0.3)'   

test_newick = Newick(test + ':0.0')
tree = BinaryTree()
tree.root = test_newick.divide()
print(tree.root.right.length)



base_length = 0.1
base_matrix=[[0.9,0.05,0.025,0.025],[0.05,0.9,0.025,0.025],[0.025,0.025,0.9,0.5],[0.025,0.025,0.5,0.9]]#the substitution matrix for length = 0.1
tree.root.left.left.value = [1,0,0,0]
tree.root.left.right.value = [1,0,0,0]
tree.root.right.value = [1,0,0,0]
print(tree.root.left.length)
pruning(base_matrix,tree.root)
