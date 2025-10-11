"""
Test runner for phylo module - Python version using biotite
"""

import numpy as np
import time
import biotite.sequence.phylo as phylo


def test_distances():
    """Test distance calculations in a tree."""
    distances = np.loadtxt("../tests/sequence/data/distances.txt", dtype=int)
    
    tree = phylo.upgma(distances)
    
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert abs(leaf.distance_to(tree.root) - dist) < 1e-6
    
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10


def test_upgma():
    """Test UPGMA algorithm."""
    # CHANGED THIS LINE - added ../
    distances = np.loadtxt("../tests/sequence/data/distances.txt", dtype=int)
    
    tree = phylo.upgma(distances)
    
    assert len(tree) == 20
    
    for i, leaf in enumerate(tree.leaves):
        assert leaf.index == i
    
    for i in range(len(tree)):
        for j in range(len(tree)):
            dist_ij = tree.get_distance(i, j)
            dist_ji = tree.get_distance(j, i)
            assert abs(dist_ij - dist_ji) < 1e-3


def test_neighbor_joining():
    """Test neighbor joining algorithm."""
    dist = np.array([
        [ 0,  5,  4,  7,  6,  8],
        [ 5,  0,  7, 10,  9, 11],
        [ 4,  7,  0,  7,  6,  8],
        [ 7, 10,  7,  0,  5,  9],
        [ 6,  9,  6,  5,  0,  8],
        [ 8, 11,  8,  9,  8,  0],
    ], dtype=np.float64)

    ref_tree = phylo.Tree(
        phylo.TreeNode(
            [
                phylo.TreeNode(
                    [
                        phylo.TreeNode(
                            [
                                phylo.TreeNode(index=0),
                                phylo.TreeNode(index=1),
                            ],
                            [1, 4],
                        ),
                        phylo.TreeNode(index=2),
                    ],
                    [1, 2],
                ),
                phylo.TreeNode(
                    [
                        phylo.TreeNode(index=3),
                        phylo.TreeNode(index=4),
                    ],
                    [3, 2],
                ),
                phylo.TreeNode(index=5),
            ],
            [1, 1, 5],
        )
    )

    test_tree = phylo.neighbor_joining(dist)

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
    
    print(f"python      {elapsed_ms}ms")