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

class Newick:
	def __init__(self,newick):
		self.newick = newick
	def divide(self):
		if ',' not in self.newick:
			length = float(self.newick.split(':')[1])
			tmp = Node(length)
			print length
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

####### use the simple test Newick string to check result
test = '((Human:0.3,Chimpanzee:0.2):0.1,Gorilla:0.3)'

test_newick = Newick(test + ':0.0')
tree = BinaryTree()
tree.root = test_newick.divide()
print(tree.root.left.right.length)
