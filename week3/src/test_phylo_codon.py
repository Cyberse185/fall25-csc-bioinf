"""
Test runner for phylo module - Codon version
"""

import time
from tree import Tree, TreeNode
from upgma import upgma
from nj import neighbor_joining
import numpy as np
from python import numpy as pnp
import numpy.pybridge


def test_distances():
    """Test distance calculations in a tree."""
    # Load the real distances.txt file as shown in instructor's hints
    distances_py: np.ndarray[int,2] = pnp.loadtxt("../tests/sequence/data/distances.txt", dtype=pnp.int64)
    # Convert to float64 for calculations
    distances = np.array(distances_py, dtype=np.float64)
    
    tree = upgma(distances)
    
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert abs(leaf.distance_to(tree.root) - dist) < 1e-6
    
    # These assertions expect 20 nodes (from distances.txt)
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10


def test_upgma():
    """Test UPGMA algorithm."""
    # Load the real distances.txt file - CHANGED THIS LINE
    distances_py: np.ndarray[int,2] = pnp.loadtxt("../tests/sequence/data/distances.txt", dtype=pnp.int64)
    distances = np.array(distances_py, dtype=np.float64)
    
    tree = upgma(distances)
    
    assert len(tree) == 20  # 20 nodes from distances.txt
    
    for i, leaf in enumerate(tree.leaves):
        assert leaf.index == i
    
    for i in range(len(tree)):
        for j in range(len(tree)):
            dist_ij = tree.get_distance(i, j)
            dist_ji = tree.get_distance(j, i)
            assert abs(dist_ij - dist_ji) < 1e-6


def test_neighbor_joining():
    """Test neighbor joining algorithm."""
    dist = np.array([
        [ 0.0,  5.0,  4.0,  7.0,  6.0,  8.0],
        [ 5.0,  0.0,  7.0, 10.0,  9.0, 11.0],
        [ 4.0,  7.0,  0.0,  7.0,  6.0,  8.0],
        [ 7.0, 10.0,  7.0,  0.0,  5.0,  9.0],
        [ 6.0,  9.0,  6.0,  5.0,  0.0,  8.0],
        [ 8.0, 11.0,  8.0,  9.0,  8.0,  0.0],
    ], dtype=np.float64)

    ref_tree = Tree(
        TreeNode(
            [
                TreeNode(
                    [
                        TreeNode(
                            [
                                TreeNode(index=0),
                                TreeNode(index=1),
                            ],
                            [1.0, 4.0],
                        ),
                        TreeNode(index=2),
                    ],
                    [1.0, 2.0],
                ),
                TreeNode(
                    [
                        TreeNode(index=3),
                        TreeNode(index=4),
                    ],
                    [3.0, 2.0],
                ),
                TreeNode(index=5),
            ],
            [1.0, 1.0, 5.0],
        )
    )

    test_tree = neighbor_joining(dist)

    for i in range(len(test_tree)):
        for j in range(len(test_tree)):
            test_dist = test_tree.get_distance(i, j)
            ref_dist = ref_tree.get_distance(i, j)
            assert abs(test_dist - ref_dist) < 1e-6


if __name__ == "__main__":
    start_time = time.time()
    
    test_distances()
    test_upgma()
    test_neighbor_joining()
    
    end_time = time.time()
    elapsed_ms = max(1, round((end_time - start_time) * 1000))
    
    print(f"codon       {elapsed_ms}ms")