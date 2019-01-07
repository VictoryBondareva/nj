import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
import networkx as nx
from nj import NJ_tree
import unittest


def read_matrix(filename):
    with open(filename, 'r') as f:
        matrix = []
        names = []
        for i, line in enumerate(f):
            if i == 0:
                names = line[:-1].split(' ')
                continue
            matrix.append(list(map(int, line[:-1].split(' '))))
    return names, matrix
# print(NJ_tree().create_tree(None))

class TestNJ_tree(unittest.TestCase):
    def test_correct_answer(self):
        for i in range(6):
            n, m = read_matrix('tests/test{}.txt'.format(i))
            tree1 = DistanceTreeConstructor().nj(DistanceMatrix(n, m))
            tree2 = NJ_tree().create_tree(n, m)
            
            self.assertTrue(nx.is_isomorphic(Phylo.to_networkx(tree1).to_undirected(), Phylo.to_networkx(tree2).to_undirected()))
    def test_empty_matrix(self):
        tree = NJ_tree().create_tree(None, None)
        self.assertEqual(tree, None)

    def test_trivial(self):
        n, m = read_matrix('tests/test_trivial.txt')
        tree = NJ_tree().create_tree(n, m)
        self.assertEqual(len(tree.clade), 0)
        
    def test_bed_matrix(self):
        n, m = read_matrix('tests/test_bad_matrix.txt')
        with self.assertRaises(ValueError):
            tree = NJ_tree().create_tree(n, m)
    
    
if __name__ == '__main__':
    unittest.main(exit=False) 

