#this script is the implementation of pruning algorithm 


base_length = 0.1
base_matrix=[[0.9,0.05,0.025,0.025],[0.05,0.9,0.025,0.025],[0.025,0.025,0.9,0.5],[0.025,0.025,0.5,0.9]]#the substitution matrix for length = 0.1
background_freq = [0.25,0.25,0.25,0.25]#assume equal probabilty for A,C,G,T



#this function does 4x4 matrix multiplication
def matrix_multiply(matrixA,matrixB):
	output = [[0 for i in range(4)]for i in range(4)]
	for i in range(4):
		for j in range(4):
			output[i][j] = matrixA[i][0]*matrixB[0][j]+matrixA[i][1]*matrixB[1][j]+matrixA[i][2]*matrixB[2][j]+matrixA[i][3]*matrixB[3][j]
	
	return output



#this function take input of base_substitution matrix and the value of t, outputp(t)
def Compute_substitution(base_matrix,length):
	relative_length = length /0.1
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


#test of Compute_substitution function

test_matrix= [[0.906563,0.045855,0.023791,0.023791],[0.045855,0.906563,0.023791,0.023791],[0.023791,0.023791,0.906563,0.045855],[0.023791,0.023791,0.045855,0.906563]]

print(sum(test_matrix[0]),sum(test_matrix[1]),sum(test_matrix[2]),sum(test_matrix[3]))

new_matrix = Compute_substitution(test_matrix,0.2)

print(new_matrix)
