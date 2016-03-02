#this test.py is to parse Newick tree file and build a binary tree data structure 

from itertools import groupby

class BinaryTree:
	def __init__(self,newick):
                '''read the newick string and return the root of the tree'''
                self.newick = newick
	        self.left = None
		self.right = None
		self.value = [0.25,0.25,0.25,0.25]#the background probability
		self.length = 0# initial length to parent is 0
	def getLeftChild(self):
		return self.left
	def getRightChild(self):
		return self.right
	def setNodeLength(self,length):#this is the length betwenn this node and its parent
		self.length = length
	def setNodeValue(self,value):
		self.value = value
        def insertLeft(self,newNode):
                self.left = newNode
        def insertRight(self,newNode):
                self.right = newNode
            
        def break_tree(self):
            	tmp = self.newick
		if ',' not in tmp:
			return False
		else:
			tab = [''.join(g) for k,g in groupby(tmp)]
			number_left = len(tab[0]) - 1# number of left parenthesis
			number_right = 0
			for i in range(len(tmp)):
				if  tmp[i] == ')':
					number_right += 1
				if number_right == number_left:
					for j in range(i+1,len(tmp)):
						if tmp[j] == ',':
							break
					break
			self.left = BinaryTree(tmp[:j] + ')')
			self.right= BinaryTree('(' + tmp[j+1:])
			return True
    
#build the tree given Newick tree
def build_tree(Newick):
	root = BinaryTree(Newick)
	if not root.break_tree():#ony one node exists in this graph
		return 
	else:
		root.break_tree()
		left_sub = root.getLeftChild()
		right_sub = root.getRightChild()
		build_tree(left_sub.newick)
		build_tree(right_sub.newick)
		return


####### use the simple test Newick string to check result		
test = '((Human:0.3,Chimpanzee:0.2):0.1,Gorilla:0.3)'   
t = BinaryTree(test)
print(t.newick)
t.break_tree()   
print(t.getRightChild().newick)
build_tree(test)
