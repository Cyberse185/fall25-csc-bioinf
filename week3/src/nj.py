# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

"""
Neighbor-Joining algorithm ported to Codon.
"""

__name__ = "phylo"
__author__ = "Patrick Kunzmann"
__all__ = ["neighbor_joining"]

from tree import Tree, TreeNode
import numpy as np
from numpy.ndarray import ndarray
from typing import List


MAX_FLOAT = np.finfo(np.float64).max


def neighbor_joining(distances: ndarray) -> Tree:
    """
    neighbor_joining(distances)
    
    Perform hierarchical clustering using the
    neighbor joining algorithm.

    In contrast to UPGMA this algorithm does not assume a constant
    evolution rate. The resulting tree is considered to be unrooted.

    Parameters
    ----------
    distances : ndarray, shape=(n,n)
        Pairwise distance matrix.

    Returns
    -------
    tree : Tree
        A rooted tree. The `index` attribute in the leaf
        TreeNode objects refer to the indices of `distances`.

    Raises
    ------
    ValueError
        If the distance matrix is not symmetric
        or if any matrix entry is below 0.
    
    Notes
    -----
    The created tree is binary except for the root node, that has three
    child nodes.

    Examples
    --------
    
    >>> distances = np.array([
    ...     [0, 1, 7, 7, 9],
    ...     [1, 0, 7, 6, 8],
    ...     [7, 7, 0, 2, 4],
    ...     [7, 6, 2, 0, 3],
    ...     [9, 8, 4, 3, 0],
    ... ])
    >>> tree = neighbor_joining(distances)
    >>> print(tree.to_newick(include_distance=False))
    (3,(2,(1,0)),4);
    """
    
    # Validation
    if distances.shape[0] != distances.shape[1] or not np.allclose(distances.T, distances):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if distances.shape[0] < 4:
        raise ValueError("At least 4 nodes are required")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")
    
    # Initialize tracking arrays
    n = distances.shape[0]
    nodes: List[TreeNode] = [TreeNode(index=i) for i in range(n)]
    
    # Track which nodes have been clustered
    is_clustered = np.full(n, False, dtype=np.uint8)
    n_rem_nodes = n - int(np.count_nonzero(is_clustered))
    
    # The divergence of a 'taxon' (relative evolution rate)
    divergence = np.zeros(n, dtype=np.float64)
    
    # Matrix for storing divergence-corrected distances
    corr_distances = np.zeros((n, n), dtype=np.float64)
    
    # Working copy of distances
    distances_v = distances.astype(np.float64, copy=True)
    
    # Main clustering loop
    while True:
        # Calculate divergence for each unclustered node
        for i in range(distances_v.shape[0]):
            if is_clustered[i]:
                continue
            dist_sum = 0.0
            for k in range(distances_v.shape[0]):
                if is_clustered[k]:
                    continue
                dist_sum += distances_v[i, k]
            divergence[i] = dist_sum
        
        # Calculate corrected distance matrix
        for i in range(distances_v.shape[0]):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                corr_distances[i, j] = (
                    float(n_rem_nodes - 2) * distances_v[i, j]
                    - divergence[i] - divergence[j]
                )
        
        # Find minimum corrected distance
        dist_min = MAX_FLOAT
        i_min = -1
        j_min = -1
        for i in range(corr_distances.shape[0]):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                dist = corr_distances[i, j]
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j
        
        # Check if all nodes have been clustered
        if i_min == -1 or j_min == -1:
            # No distance found -> all leaf nodes are clustered
            break
        
        # Cluster the nodes with minimum distance
        node_dist_i = 0.5 * (
            distances_v[i_min, j_min]
            + 1.0 / float(n_rem_nodes - 2) * (divergence[i_min] - divergence[j_min])
        )
        node_dist_j = 0.5 * (
            distances_v[i_min, j_min]
            + 1.0 / float(n_rem_nodes - 2) * (divergence[j_min] - divergence[i_min])
        )
        
        if n_rem_nodes > 3:
            # Clustering is not finished
            # Create a node with two children
            nodes[i_min] = TreeNode(
                [nodes[i_min], nodes[j_min]],
                [node_dist_i, node_dist_j]
            )
            # Mark position j_min as clustered
            is_clustered[j_min] = True
        else:
            # Clustering is finished
            # Combine last three nodes into root node
            is_clustered[i_min] = True
            is_clustered[j_min] = True
            # Find the index of the remaining node (other than i_min and j_min)
            k = -1
            for idx in range(n):
                if not is_clustered[idx]:
                    k = idx
                    break
            node_dist_k = 0.5 * (
                distances_v[i_min, k] + distances_v[j_min, k]
                - distances_v[i_min, j_min]
            )
            root = TreeNode(
                [nodes[i_min], nodes[j_min], nodes[k]],
                [node_dist_i, node_dist_j, node_dist_k]
            )
            # Clustering is finished -> return tree
            return Tree(root)
        
        # Update distance matrix
        # Calculate distances of new node to all other nodes
        for k in range(distances_v.shape[0]):
            if not is_clustered[k] and k != i_min:
                dist = 0.5 * (
                    distances_v[i_min, k] + distances_v[j_min, k]
                    - distances_v[i_min, j_min]
                )
                distances_v[i_min, k] = dist
                distances_v[k, i_min] = dist
        
        # Update the count of remaining nodes
        n_rem_nodes = n - int(np.count_nonzero(is_clustered))
    
    # Should never reach here, but return a dummy tree just in case
    raise ValueError("Clustering failed unexpectedly")