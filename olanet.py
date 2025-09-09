import argparse
import csv
import sys
import random
from copy import copy
from typing import List, Dict, Tuple, Set, Optional
import multiprocessing
import os

import treeswift
from treeswift import Tree, Node


def read_tree_file(file_path: str) -> List[Tree]:
    """Reads a tree file in Newick or Nexus format and returns a list of trees."""
    try:
        trees = treeswift.read_tree_newick(file_path)
        if not isinstance(trees, list):
            return [trees]
        return trees
    except Exception as e_newick:
        try:
            nexus_data = treeswift.read_tree_nexus(file_path)
            return [val for val in nexus_data.values() if isinstance(val, Tree)]
        except Exception as e_nexus:
            print(f"Error: Could not parse tree file '{file_path}'.", file=sys.stderr)
            print(f"Newick parser error: {e_newick}", file=sys.stderr)
            print(f"Nexus parser error: {e_nexus}", file=sys.stderr)
            sys.exit(1)


def read_dates_csv(file_path: str) -> Dict[str, str]:
    """Reads a CSV file with taxon and date columns and returns a dictionary."""
    dates_map: Dict[str, str] = {}
    try:
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            _ = next(reader, None)  # Skip header
            for row in reader:
                if len(row) == 2:
                    taxon, date = row
                    dates_map[taxon.strip()] = date.strip()
                else:
                    print(f"Warning: Skipping invalid row in '{file_path}': {row}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Dates file not found at '{file_path}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Could not read or parse dates file '{file_path}': {e}", file=sys.stderr)
        sys.exit(1)
    return dates_map


def prune_trees_to_common_taxa(trees: List[Tree]) -> List[Tree]:
    """
    Finds the common set of taxa among all trees, prunes each tree to that common set, and returns pruned trees.
    """
    # Start with the taxa from the first tree
    common_taxa = set(trees[0].labels())

    # Find the intersection of taxa with all other trees
    for tree in trees[1:]:
        common_taxa.intersection_update(set(tree.labels()))

    if not common_taxa:
        print("Error: No common taxa found across trees found - exiting.")
        sys.exit(1)
    if len(common_taxa) < 3:
        print("Error: Very few common taxa across trees (less than 3) - exiting.")
        sys.exit(1)

    print(f"Found {len(common_taxa)} common taxa across all trees.")

    # Prune each tree to the common set of taxa
    pruned_trees = []
    for i, tree in enumerate(trees):
        pruned_tree = tree.extract_tree_with(common_taxa)

        pruned_trees.append(pruned_tree)
    return pruned_trees


def check_taxa_consistency(trees: List[Tree]):
    """Checks if all trees have the same set of leaf labels."""
    if not trees:
        return
    first_tree_labels = set(trees[0].labels())
    for i in range(1, len(trees)):
        if set(trees[i].labels()) != first_tree_labels:
            print(f"Error: Trees have inconsistent leaf labels. Tree {i+1} does not match the first tree.", file=sys.stderr)
            sys.exit(1)
    print("All trees have a consistent set of leaf labels.")


def hamming_ola_distance(vectors: List[List[int]]) -> int:
    """
    Computes the number of indices where at least two OLA vectors disagree.
    """
    if not vectors or len(vectors) < 2:
        return 0

    # Check that all vectors have the same length
    vector_iter = iter(vectors)
    first_len = len(next(vector_iter))
    if not all(len(vec) == first_len for vec in vector_iter):
        print("Error: OLA vectors must have the same length for comparison.", file=sys.stderr)
        return -1

    num_indices = first_len
    disagreement_count = 0
    for i in range(num_indices):
        # Get all values at the current index
        values_at_index = {v[i] for v in vectors}
        # If the set of values has more than one element, there's a disagreement
        if len(values_at_index) > 1:
            disagreement_count += 1

    return disagreement_count


def corrected_ola_distance(vectors: List[List[int]]) -> int:
    """
    Computes the "corrected" number of mismatched indices (reticulation number).
    An index i is added to a set of mismatches M if
        (i) there is disagreement across vectors, or
        (ii) all vectors store -j, where j is in M.
    """
    if not vectors or len(vectors) < 2:
        return 0

    first_len = len(vectors[0])
    if not all(len(vec) == first_len for vec in vectors):
        print("Error: OLA vectors must have the same length for comparison.", file=sys.stderr)
        return 1

    mismatches = [0] * (first_len + 1)
    count = 0
    for i in range(first_len):
        leaf_ind = i + 1
        values_at_index = {v[i] for v in vectors}
        # If the set of values has more than one element, there's a disagreement
        if len(values_at_index) > 1:
            count += 1
            mismatches[leaf_ind] = 1
        else:
            consensus_value = values_at_index.pop()
            if consensus_value < 0 and mismatches[-consensus_value] > 0:
                # If the consensus match is at a parent of a mismatched index.
                count += 1
                mismatches[leaf_ind] = 1
    return count


def ola_from_binary_tree(tree: Tree, leaf_ordering: List[str]) -> List[int]:
    """
    Converts a tree to an OLA vector given a leaf ordering.
    Implementation based on "Vector encoding of phylogenetic trees by ordered leaf attachment"
    by Richman, Zhang, and Matsen (2025).
    """
    n = len(leaf_ordering)
    if n <= 1:
        return []

    # Work on a copy of the tree to avoid modifying the original
    tree_copy = copy(tree)
    label_to_int = {label: i for i, label in enumerate(leaf_ordering)}

    # --- Phase 1: Assign Canonical Internal Node Labels (Algorithm 1) ---
    # We attach the labels directly to the nodes of the copied tree.
    for node in tree_copy.traverse_postorder():
        if node.is_leaf():
            node.cf_label = label_to_int[node.label]
        else:
            node.cf_label = min(child.cf_label for child in node.children)

    for node in tree_copy.traverse_preorder():
        if node.is_leaf():
            node.canonical_label = label_to_int[node.label]
        else:
            cs_label = max(child.cf_label for child in node.children)
            node.canonical_label = -cs_label

    # --- Phase 2: Generate OLA Vector (Algorithm 2) ---
    leaf_nodes = {label_to_int[n.label]: n for n in tree_copy.traverse_leaves()}
    ola_vector = [0] * (n - 1)

    # Deconstruct the tree by finding the sister of each leaf in reverse order
    for i in range(n - 1, 0, -1):
        leaf_node = leaf_nodes[i]
        parent = leaf_node.parent

        sister_node = next(child for child in parent.children if child != leaf_node)
        ola_vector[i - 1] = sister_node.canonical_label

        # Modify the tree structure for the next iteration
        grandparent = parent.parent
        if grandparent:
            grandparent.remove_child(parent)
            grandparent.add_child(sister_node)
        else:
            # The parent was the root, so the sister becomes the new root
            tree_copy.root = sister_node
            sister_node.parent = None

    return ola_vector


def preprocess_tree(tree: Tree, leaf_ordering: List[str]) -> Tuple[List[int], List[bool]]:
    """
    Non-binary tree preprocessing: generate an initial OLA vector and identify multifurcations for a single tree.
    """
    tree.suppress_unifurcations()  # We do not allow unifurcations.
    n = len(leaf_ordering)
    ola_vector = [0] * (n - 1)
    multi_vector = [False] * (n - 1)

    tree_copy = copy(tree)
    label_to_int = {label: i for i, label in enumerate(leaf_ordering)}

    for node in tree_copy.traverse_postorder():
        if node.is_leaf():
            node.min_leaf_idx = label_to_int[node.label]
            node.canonical_idx = label_to_int[node.label]
        else:
            children_min_indices = sorted([c.min_leaf_idx for c in node.children])
            node.min_leaf_idx = children_min_indices[0]
            node.canonical_idx = -children_min_indices[1]

    leaf_nodes = {label_to_int[leaf.label]: leaf for leaf in tree_copy.traverse_leaves()}

    for i in range(n - 1, 0, -1):
        leaf_node = leaf_nodes[i]
        parent = leaf_node.parent
        is_multi = len(parent.children) > 2
        multi_vector[i - 1] = is_multi

        if is_multi:
            ola_vector[i - 1] = parent.canonical_idx
            parent.remove_child(leaf_node)
        else:
            sibling = next(c for c in parent.children if c != leaf_node)
            ola_vector[i - 1] = sibling.canonical_idx
            
            grandparent = parent.parent
            if grandparent:
                grandparent.remove_child(parent)
                grandparent.add_child(sibling)
            else:
                tree_copy.root = sibling
                sibling.parent = None
    
    return ola_vector, multi_vector


def _insert_leaf(place_node: Node, leaf_label: str, tree: Tree) -> Tuple[Node, Node]:
    """Inserts a new leaf above a given place_node, returning (leaf, new_internal_node)"""
    new_leaf = treeswift.Node(label=leaf_label)
    new_internal_node = treeswift.Node()
    
    grandparent = place_node.parent
    if grandparent:
        grandparent.remove_child(place_node)
        grandparent.add_child(new_internal_node)
    else:
        tree.root = new_internal_node
    
    new_internal_node.add_child(place_node)
    new_internal_node.add_child(new_leaf)
    return new_leaf, new_internal_node


def _add_leaf_to_resolved_tree(leaf_id: int, placement_id: int, tree_id: int, resolved_trees: List[Tree],
                               node_maps: List[Dict[int, Node]], m_sets: List[List[Set[int]]], m_roots: List[List[int]],
                               leaf_of: List[Dict[int, int]], leaf_label: str):
    place_node = node_maps[tree_id][placement_id]
    new_leaf, new_internal = _insert_leaf(place_node, leaf_label, resolved_trees[tree_id])
    
    node_maps[tree_id][leaf_id] = new_leaf
    node_maps[tree_id][-leaf_id] = new_internal
    m_sets[tree_id][leaf_id] = {leaf_id, -leaf_id, placement_id}
    m_roots[tree_id][leaf_id] = -leaf_id
    
    y = leaf_of[tree_id].get(placement_id)
    if y and y > 0:
        m_sets[tree_id][y].discard(placement_id)
        m_sets[tree_id][y].add(-leaf_id)
        leaf_of[tree_id][-leaf_id] = y
    leaf_of[tree_id][leaf_id] = leaf_of[tree_id][placement_id] = leaf_id


def _determine_consensus_placement(leaf_id: int, placements: Dict[int, int], unresolved_indices: List[int],
                                   ola_vectors: List[List[int]], m_sets: List[List[Set[int]]],
                                   mismatch_neg: Set[int]) -> Optional[int]:
    if len(set(placements.values())) > 1:
        mismatch_neg.add(-leaf_id)
        return None
    if placements:
        return list(placements.values())[0]
    
    intersection = None
    for tree_id in unresolved_indices:
        multi_id = ola_vectors[tree_id][leaf_id - 1]
        current_set = m_sets[tree_id][abs(multi_id)]
        if intersection is None:
            intersection = current_set.copy()
        else:
            intersection.intersection_update(current_set)
    if intersection:
        intersection.difference_update(mismatch_neg)

    if not intersection:
        mismatch_neg.add(-leaf_id)
        return None
    else:
        return max(intersection, key=abs)


def _add_leaf_to_unresolved_tree(leaf_id: int, consensus_placement: Optional[int], tree_id: int,
                                 resolved_trees: List[Tree], node_maps: List[Dict[int, Node]],
                                 ola_vectors: List[List[int]], m_sets: List[List[Set[int]]], m_roots: List[List[int]],
                                 leaf_of: List[Dict[int, int]], leaf_label: str):
    mult_id = abs(ola_vectors[tree_id][leaf_id - 1])
    placement = consensus_placement if consensus_placement is not None and consensus_placement in m_sets[tree_id][mult_id] else m_roots[tree_id][mult_id]
    
    place_node = node_maps[tree_id][placement]
    new_leaf, new_internal = _insert_leaf(place_node, leaf_label, resolved_trees[tree_id])

    node_maps[tree_id][leaf_id] = new_leaf
    node_maps[tree_id][-leaf_id] = new_internal
    
    m_sets[tree_id][mult_id].update({leaf_id, -leaf_id})
    leaf_of[tree_id][leaf_id] = mult_id
    
    if placement == m_roots[tree_id][mult_id]:
        m_roots[tree_id][mult_id] = -leaf_id
        
    y = leaf_of[tree_id].get(placement)
    if y and y != mult_id:
        m_sets[tree_id][y].discard(placement)
        m_sets[tree_id][y].add(-leaf_id)
        leaf_of[tree_id][-leaf_id] = y
        leaf_of[tree_id].pop(placement, None)


def resolve_multifurcations(trees: List[Tree], leaf_ordering: List[str]) -> Tuple[List[Tree], List[List[int]]]:
    """
    Optimally (jointly) resolves trees given a fixed leaf ordering.
    """
    n = len(leaf_ordering)
    num_trees = len(trees)
    int_to_label = {i: label for i, label in enumerate(leaf_ordering)}

    resolved_trees = []
    for _ in trees:
        t_prime = treeswift.Tree()
        t_prime.root.set_label(int_to_label[0])
        resolved_trees.append(t_prime)

    ola_vectors = [[] for _ in trees]
    multi_vectors = [[] for _ in trees]
    for i, t in enumerate(trees):
        ola_vectors[i], multi_vectors[i] = preprocess_tree(t, leaf_ordering)

    mismatch_neg: Set[int] = set()
    m_sets = [[set() for _ in range(n)] for _ in range(num_trees)]
    m_roots = [[0 for _ in range(n)] for _ in range(num_trees)]
    leaf_of = [{} for _ in range(num_trees)]
    node_maps = []
    for t_prime in resolved_trees:
        node_map = {0: t_prime.root}
        node_maps.append(node_map)

    for i in range(1, n):
        leaf_label = int_to_label[i]
        
        placements = {}
        unresolved_indices = []
        for tree_id in range(num_trees):
            if multi_vectors[tree_id][i - 1]:
                unresolved_indices.append(tree_id)
            else:
                mult_placement = ola_vectors[tree_id][i - 1]
                if mult_placement < 0:
                    placements[tree_id] = m_roots[tree_id][-mult_placement]
                else:
                    placements[tree_id] = mult_placement

        for tree_id, p_idx in placements.items():
            _add_leaf_to_resolved_tree(i, p_idx, tree_id, resolved_trees, node_maps, m_sets, m_roots, leaf_of, leaf_label)

        p = _determine_consensus_placement(i, placements, unresolved_indices, ola_vectors, m_sets, mismatch_neg)

        for tree_id in unresolved_indices:
            _add_leaf_to_unresolved_tree(i, p, tree_id, resolved_trees, node_maps, ola_vectors, m_sets, m_roots, leaf_of, leaf_label)

    final_ola_vectors = [ola_from_binary_tree(t, leaf_ordering) for t in resolved_trees]
    return resolved_trees, final_ola_vectors


def _process_permutation(trees: List[Tree], ordering: List[str]) -> int:
    """Helper for multiprocessing."""
    _, resolved_vectors = resolve_multifurcations(trees, ordering)
    return corrected_ola_distance(resolved_vectors)


def main():
    parser = argparse.ArgumentParser(
        description='OLA-Net: Optimally resolve multifurcations in a set of trees and compute the reticulation number.'
    )

    parser.add_argument(
        '-t', '--trees',
        required=True,
        help='A file with phylogenetic trees (Newick or Nexus format).'
    )
    parser.add_argument(
        '-b', '--support',
        type=float,
        default=70,
        help='Collapse branches with support below the threshold (default: 70).'
    )
    parser.add_argument(
        '-s', '--short',
        type=float,
        default=1e-9,
        help='A threshold for collapsing short branches (default: 1e-9).'
    )
    parser.add_argument(
        '--dates',
        help='An optional CSV file of type "taxon,date" to define leaf ordering. Dates must be in YYYY-mm-dd format.'
    )
    parser.add_argument(
        '-p', '--permutations',
        type=int,
        default=1000,
        help='Number of random permutations to try if --dates is not provided (default: 1000).'
    )
    parser.add_argument(
        '--cores',
        type=int,
        default=os.cpu_count(),
        help='Number of CPU cores to use for parallel processing (default: all available cores).'
    )

    args = parser.parse_args()

    print("Reading tree file...")
    trees = read_tree_file(args.trees)
    print(f"{len(trees)} tree(s) loaded successfully from '{args.trees}'.")
    if len(trees) < 2:
        parser.error("Error: Less than two trees in the trees file - exiting.")

    trees = prune_trees_to_common_taxa(trees)
    all_leaf_labels = list(trees[0].labels(internal=False))

    # Collapse short and low-supported branches
    print(f'Collapsing branches with support below {args.support} and length below {args.short} '
          f'(if branch lengths and/or support values are specified in trees)')
    for i, tree in enumerate(trees):
        print(f'Tree {i+1}')
        print('\tInitial # of internal nodes:', tree.num_nodes(leaves=False, internal=True))
        if tree.edge_length_sum() > 0:
            # Collapse short branches if the tree has edge lengths.
            tree.collapse_short_branches(args.short)
        print('\tAfter short branches collapsed:', tree.num_nodes(leaves=False, internal=True))
        tree.contract_low_support(args.support, terminal=False, internal=True)
        print('\tAfter low-support collapsed:', tree.num_nodes(leaves=False, internal=True))

    if args.dates:
        print("Reading dates file for leaf ordering...")
        dates_map = read_dates_csv(args.dates)
        
        if not set(dates_map.keys()).issuperset(set(all_leaf_labels)):
            print(set(dates_map.keys()))
            print(all_leaf_labels)
            print("Error: Taxa in dates file do not match taxa in trees.", file=sys.stderr)
            sys.exit(1)

        sorted_taxa = sorted(dates_map.items(), key=lambda item: item[1])
        date_leaf_ordering = [taxon for taxon, date in sorted_taxa]
        
        print(f"Using the date-defined leaf ordering for multifurcation resolution...")
        _, resolved_vectors = resolve_multifurcations(trees, date_leaf_ordering)

        distance = corrected_ola_distance(resolved_vectors)
        print(f"\nEstimated reticulation number: {distance}")
        # ham_dist = hamming_ola_distance(resolved_vectors)
        # print(f"\nHamming OLA distance: {ham_dist}")

    else:
        print(f"Generating {args.permutations} random leaf orderings...")
        leaf_orderings = []
        base_labels = all_leaf_labels
        for _ in range(args.permutations):
            random.shuffle(base_labels)
            leaf_orderings.append(list(base_labels))
        
        print(f"Processing {len(leaf_orderings)} permutations using {args.cores} cores...")
        
        pool_args = [(trees, o) for o in leaf_orderings]
        
        with multiprocessing.Pool(processes=args.cores) as pool:
            results = pool.starmap(_process_permutation, pool_args)

        min_distance = min(results) if results else float('inf')
        print(f"\nMinimum Estimated Reticulation Number: {min_distance}")


if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    main()
