from __future__ import division
import re

whatever = '((Human:0.3,Chimpanzee:0.2):0.1,Gorilla:0.3)'

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

print parseNewick(whatever, name='Root')
