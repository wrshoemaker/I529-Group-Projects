#!/usr/bin/env python
from __future__ import division
import itertools, collections, os, math, re, argparse, operator
from collections import defaultdict, Iterable
import collections

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'
#treeIN = mydir + '/' + 'tree_file.txt'
#alignIN = mydir + '/' + 'alignment_file.txt'


def readNewick(fileNewick):
    listNewick = []
    with open(fileNewick, "r") as f:
        line = f.read().strip()
        listNewick.append(line)
    return listNewick

def readAlignment(fileAlignment):
    '''
    Takes the file containing the sequences and returns a
    nested list containing lists of ['name', 'sequence']
    '''
    listAlignment = []
    with open(fileAlignment, "r") as f:
        for line in f:
            line = line.split(':')
            line = map(str.strip, line)
            listAlignment.append(line)
    return listAlignment

def formatAlignment(listAlignment):
    '''
    Takes the output from readAlignment and returns a tuple containing to lists;
    the first list contains sequence names, the second list contains tuples
    of the aligned nucleotides. ex (['seq1','seq2'], [('A','C'), ('A','G')])
    Uninformative sites are removed.
    '''
    # get just the sequences
    sequences = [ x[1].upper() for x in listAlignment ]
    names = [ x[0] for x in listAlignment ]
#    print sequences
    alignedSeqs = zip(*sequences)
#    print alignedSeqs
#    for x in alignedSeqs:
#        if len(set(x)) < 2:
#           alignedSeqs.remove(x)
    return (names, alignedSeqs)

def parseNewick(newick, name='Root'):
    tree = []
    i = 0
    while i < len(newick) - 1:
        i = i + 1
        if startsClade(newick[i]):
            clade = getClade(newick[i:])
            node_name = getWord(newick[i + len(clade)])
            tree.append(parseNewick(clade, node_name))
            i = i + len(clade) + len(node_name)

        elif startsLeaf(newick[i]):
            leaf = getWord(newick[i:])

            if newick.count(')') > 1:
                test_what = newick.rsplit(')', 2)
                str_list = filter(None, test_what)
                return_values = getWord_test(str_list[-1])
            else:
                return_values = getWord_test(newick[i:])
            tree.append(return_values)
            i = i + len(leaf)
    return tree

def startsLeaf(c):
    return c.isalpha()

def startsClade(c):
    return c == '('

def endsClade(c):
    return c == ')'

def getClade(newick):
    desc = 0
    for i, c in enumerate(newick):
        if startsClade(c):
            desc += 1
        elif endsClade(c):
            desc -= 1
        if desc == 0:
            return newick[:i+1]

def getWord(string):
    for i, c in enumerate(string):
        if not c.isalnum():
            return string[:i]
    return string

def getWord_test(string):
    stringSplit = re.split(r'[(),:]+', string)
    stringSplit = filter(None, stringSplit)
    if len(stringSplit) == 4:
        return stringSplit[:-2]
    else:
        return stringSplit

###### prunning algorithm
## row1, row2, etc
###[[A, C, G, T]]
## unit of evolutionary time


def matrix_multiply(a,b):
    zip_b = zip(*b)
    return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b))
             for col_b in zip_b] for row_a in a]


#this function take input of base_substitution matrix and the value of t, outputp(t)
def Compute_substitution(base_matrix, length):
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

def getProb(fromBase, toBase, length):
    toBase = toBase.upper()
    fromBase = fromBase.upper()
    refDict = {'A':0, 'C':1, 'G':2, 'T':3}
    toIndex = refDict[toBase]
    fromIndex = refDict[fromBase]
    matrix = Compute_substitution(length)
    return matrix[fromIndex][toIndex]


def prunningAlgorithm(node1,node2,matrix1,matrix2):
    # node1, node2 are two children of one node
    node1_A, node1_C, node1_G, node1_T = node1[0], node1[1], node1[2], node1[3]
    node2_A, node2_C, node2_G, node2_T = node2[0], node2[1], node2[2], node2[3]
    #print node2_A
    # freq of parental node is A
    nodeP_A = (node1_A * matrix1[0][0] + node1_C * matrix1[0][1] + \
        node1_G * matrix1[0][2] + node1_T * matrix1[0][3]) * \
        (node2_A * matrix2[0][0] + node2_C * matrix2[0][1] + \
        node2_G * matrix2[0][2] + node2_T * matrix2[0][3])
    # freq of parental node is C
    nodeP_C = (node1_A * matrix1[1][0] + node1_C * matrix1[1][1] + \
        node1_G * matrix1[1][2] + node1_T * matrix1[1][3]) * \
        (node2_A * matrix2[1][0] + node2_C * matrix2[1][1] + \
        node2_G * matrix2[1][2] + node2_T * matrix2[1][3])
    # freq of parental node is G
    nodeP_G = (node1_A * matrix1[2][0] + node1_C * matrix1[2][1] + \
        node1_G * matrix1[2][2] + node1_T * matrix1[2][3]) * \
        (node2_A * matrix2[2][0] + node2_C * matrix2[2][1]+\
        node2_G * matrix2[2][2] + node2_T * matrix2[2][3])
    # freq of parental node is T
    nodeP_T = (node1_A * matrix1[3][0] + node1_C * matrix1[3][1] + \
        node1_G * matrix1[3][2] + node1_T * matrix1[3][3]) * \
        (node2_A * matrix2[3][0] + node2_C * matrix2[3][1] + \
        node2_G * matrix2[3][2] + node2_T * matrix2[3][3])
    nodeP = [nodeP_A,nodeP_C,nodeP_G,nodeP_T]
    return nodeP



def flatten(coll):
    for i in coll:
        if isinstance(i, Iterable) and not isinstance(i, basestring):
            for subc in flatten(i):
                yield subc
        else:
            yield i

def siteToDummy(site):
    stateList = []
    for sample in site:
        if sample == 'A':
            stateList.append([1,0,0,0])
        elif sample == 'C':
            stateList.append([0,1,0,0])
        elif sample == 'G':
            stateList.append([0,0,1,0])
        else:
            stateList.append([0,0,0,1])
    return stateList

def recursivePrunning(newickTree, alignedSeqs, subMatrix):
    '''
    This function takes the newick tree, the aligned sequences, and the substitution
    matrix as arguments. It returns a nested list, with each element containing the
    likelihood of observing a given nucleotide (list order = A, C, G, T) for
    the alignment given the evolutionary distance supplied by the phylogeny
    '''
    numberTaxa = len(alignedSeqs[0])
    # go through all n-1 comparisons (n = # nodes)
    flattenedNewick = list(flatten(newickTree))
    probSiteGivenTree = []
    for seq in alignedSeqs[1]:
        seqDummy = siteToDummy(seq)
        currentPrunning = []
        for x in range (0, numberTaxa-1):
            if x == 0:
                matrix1 = Compute_substitution(subMatrix,float(flattenedNewick[1]))
                matrix2 = Compute_substitution(subMatrix,float(flattenedNewick[3]))
                #print siteToDummy(alignedSeqs[1])
                currentPrunning = prunningAlgorithm(seqDummy[0],seqDummy[1], matrix1, matrix2)
            else:
                # matrix 1 is from the previously merged nodes
                matrix1 = Compute_substitution(subMatrix,float(flattenedNewick[ x + ( 3 * (numberTaxa-2) ) ]))
                # matrix 2 is from the outer node
                matrix2 = Compute_substitution(subMatrix,float(flattenedNewick[ x + ( 5 * (numberTaxa-2) ) ]))
                node1 = currentPrunning
                node2 = seqDummy[x + 1]
                currentPrunning = prunningAlgorithm(node1,node2, matrix1, matrix2)
        probSiteGivenTree.append(currentPrunning)
    return probSiteGivenTree


subMatrix = [[0.9,0.05,0.025,0.025],[0.05,0.9,0.025,0.025], \
        [0.025,0.025,0.9,0.05],[0.025,0.025,0.05,0.9]]

parser = argparse.ArgumentParser()
parser.add_argument('-t','--tree_file',required=True)
parser.add_argument('-a','--alignment_file',required=True)
args = parser.parse_args()
treeIN = args.tree_file
alignIN = args.alignment_file

# test seqs
test_read = readAlignment(alignIN)

alignedSequences =  formatAlignment(test_read)


#	print(right_matrixi test tree
newickTest = readNewick(treeIN)[0]
test_newick = parseNewick(newickTest, name='Root')


prob_set = recursivePrunning(test_newick, alignedSequences, subMatrix)
prob = 1
for i in prob_set:
	prob *= sum(i)
print float(prob)
