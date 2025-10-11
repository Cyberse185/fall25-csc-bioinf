# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

"""
UPGMA algorithm ported to Codon.
"""

__name__ = "phylo"
__author__ = "Patrick Kunzmann"
__all__ = ["upgma"]

from tree import Tree, TreeNode
import numpy as np
from numpy.ndarray import ndarray
from typing import List


MAX_FLOAT = np.finfo(np.float64).max


def upgma(distances: ndarray) -> Tree:
    """
    upgma(distances)
    
    Perform hierarchical clustering using the
    unweighted pair group method with arithmetic mean (UPGMA).
    
    This algorithm produces leaf nodes with the same distance to the
    root node. In the context of evolution this means a constant 
    evolution rate (molecular clock).

    Parameters
    ----------
    distances : ndarray, shape=(n,n)
        Pairwise distance matrix.

    Returns
    -------
    tree : Tree
        A rooted binary tree. The `index` attribute in the leaf
        TreeNode objects refer to the indices of `distances`.

    Raises
    ------
    ValueError
        If the distance matrix is not symmetric
        or if any matrix entry is below 0.

    Examples
    --------
    
    >>> distances = np.array([
    ...     [0, 1, 7, 7, 9],
    ...     [1, 0, 7, 6, 8],
    ...     [7, 7, 0, 2, 4],
    ...     [7, 6, 2, 0, 3],
    ...     [9, 8, 4, 3, 0],
    ... ])
    >>> tree = upgma(distances)
    >>> print(tree.to_newick(include_distance=False))
    ((4,(3,2)),(1,0));
    """
    
    # Validation
    if distances.shape[0] != distances.shape[1] or not np.allclose(distances.T, distances):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")
    
    # Initialize tracking arrays
    n = distances.shape[0]
    nodes: List[TreeNode] = [TreeNode(index=i) for i in range(n)]
    
    # Track which nodes have been clustered
    is_clustered = np.full(n, False, dtype=np.uint8)
    
    # Number of indices in current node (for proportional averaging)
    # Changed to float64 to avoid type casting issues
    cluster_size = np.ones(n, dtype=np.float64)
    
    # Distance of each node from leaf nodes
    node_heights = np.zeros(n, dtype=np.float64)
    
    # Working copy of distances
    distances_v = distances.astype(np.float64, copy=True)
    
    # Track the root index
    root_idx = 0
    
    # Main clustering loop
    while True:
        # Find minimum distance
        dist_min = MAX_FLOAT
        i_min = -1
        j_min = -1
        
        for i in range(distances_v.shape[0]):
            if is_clustered[i]:
                continue
            for j in range(i):
                if is_clustered[j]:
                    continue
                dist = distances_v[i, j]
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j
        
        if i_min == -1 or j_min == -1:
            # No distance found -> all leaf nodes are clustered
            break
        
        # Cluster the nodes with minimum distance
        height = dist_min / 2.0
        nodes[i_min] = TreeNode(
            [nodes[i_min], nodes[j_min]],
            [height - node_heights[i_min], height - node_heights[j_min]]
        )
        node_heights[i_min] = height
        root_idx = i_min
        
        # Mark position j_min as clustered
        is_clustered[j_min] = True
        
        # Calculate arithmetic mean distances of child nodes
        # as distances for new node and update matrix
        for k in range(distances_v.shape[0]):
            if not is_clustered[k] and k != i_min:
                mean = (
                    (distances_v[i_min, k] * cluster_size[i_min] + 
                     distances_v[j_min, k] * cluster_size[j_min])
                    / (cluster_size[i_min] + cluster_size[j_min])
                )
                distances_v[i_min, k] = mean
                distances_v[k, i_min] = mean
        
        # Update cluster size of new node
        cluster_size[i_min] = cluster_size[i_min] + cluster_size[j_min]
    
    # The root node is at root_idx
    return Tree(nodes[root_idx])