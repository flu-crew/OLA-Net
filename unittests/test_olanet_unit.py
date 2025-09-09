import unittest
import treeswift
from olanet import resolve_multifurcations, ola_from_binary_tree, hamming_ola_distance, corrected_ola_distance


class TestOLAVectors(unittest.TestCase):

    def test_ola_vector_1(self):
        tree = treeswift.read_tree_newick("(0,(1,(2,(3,4))));")
        leaf_ordering = ['0', '1', '2', '3', '4']
        expected_ola_vector = [0, 1, 2, 3]
        self.assertEqual(ola_from_binary_tree(tree, leaf_ordering), expected_ola_vector)

    def test_ola_vector_2(self):
        tree = treeswift.read_tree_newick("((0,((1,3),4)),2);")
        leaf_ordering = ['0', '1', '2', '3', '4']
        expected_ola_vector = [0, -1, 1, -3]
        self.assertEqual(ola_from_binary_tree(tree, leaf_ordering), expected_ola_vector)

    def test_ola_vector_3(self):
        tree = treeswift.read_tree_newick("((A,B),(C,D));")
        leaf_ordering = ['C', 'A', 'D', 'B']
        expected_ola_vector = [0, 0, 1]
        self.assertEqual(ola_from_binary_tree(tree, leaf_ordering), expected_ola_vector)

    def test_ola_distance(self):
        tree1 = treeswift.read_tree_newick("(((0,2),3),1);")
        tree2 = treeswift.read_tree_newick("(((1,2),3),0);")
        leaf_ordering = ['0', '1', '2', '3']
        expected_ola_vector1 = [0, 0, -2]
        expected_ola_vector2 = [0, 1, -2]
        expected_ham_dist = 1
        expected_cor_dist = 2
        ola_vector1 = ola_from_binary_tree(tree1, leaf_ordering)
        ola_vector2 = ola_from_binary_tree(tree2, leaf_ordering)
        self.assertEqual(ola_vector1, expected_ola_vector1)
        self.assertEqual(ola_vector2, expected_ola_vector2)
        self.assertEqual(hamming_ola_distance([ola_vector1, ola_vector2]), expected_ham_dist)
        self.assertEqual(corrected_ola_distance([ola_vector1, ola_vector2]), expected_cor_dist)

    def test_ola_multifurcated1(self):
        tree1 = treeswift.read_tree_newick("(0,(1,2,3,4));")
        tree2 = treeswift.read_tree_newick("(1,2,3,(0,4));")
        leaf_ordering = ['0', '1', '2', '3', '4']
        expected_ola_vector1 = [0, 1, 2, -2]
        expected_ola_vector2 = [0, 1, 2, 0]
        (t1_res, t2_res), (ola1, ola2) = resolve_multifurcations([tree1, tree2], leaf_ordering)
        self.assertEqual(ola_from_binary_tree(t1_res, leaf_ordering), expected_ola_vector1)
        self.assertEqual(ola1, expected_ola_vector1)
        self.assertEqual(ola_from_binary_tree(t2_res, leaf_ordering), expected_ola_vector2)
        self.assertEqual(ola2, expected_ola_vector2)
        self.assertEqual(corrected_ola_distance([ola1, ola2]), 1)

    def test_ola_multifurcated2(self):
        tree1 = treeswift.read_tree_newick("((1,2,3),(4,5,6),(7,8,9));")
        tree2 = treeswift.read_tree_newick("((((1,2),(3,4,5,6)),(7,8)),9);")
        leaf_ordering = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
        expected_ola_vector1 = [0, -1, -2, 3, 4, -2, 6, -7]
        expected_ola_vector2 = [0, -1,  2, 3, 4, -2, 6, -6]
        (t1_res, t2_res), (ola1, ola2) = resolve_multifurcations([tree1, tree2], leaf_ordering)
        self.assertEqual(ola_from_binary_tree(t1_res, leaf_ordering), expected_ola_vector1)
        self.assertEqual(ola1, expected_ola_vector1)
        self.assertEqual(ola_from_binary_tree(t2_res, leaf_ordering), expected_ola_vector2)
        self.assertEqual(ola2, expected_ola_vector2)
        self.assertEqual(corrected_ola_distance([ola1, ola2]), 2)


if __name__ == '__main__':
    unittest.main()
