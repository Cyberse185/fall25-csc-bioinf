# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

"""
Phylogenetic tree data structures ported to Codon.
"""

__name__ = "phylo"
__author__ = "Patrick Kunzmann, Tom David Müller"
__all__ = ["Tree", "TreeNode", "as_binary"]

from numpy.ndarray import ndarray
import numpy as np
from typing import Optional, List, Tuple


# Add set hash extension for Codon
@extend
class set:
    def __hash__(self):
        MAX = int.MAX
        MASK = 2 * MAX + 1
        n = len(self)
        h = 1927868237 * (n + 1)
        h &= MASK
        for x in self:
            hx = hash(x)
            h ^= (hx ^ (hx << 16) ^ 89869747) * 3644798167
            h &= MASK
        h = h * 69069 + 907133923
        h &= MASK
        if h > MAX:
            h -= MASK + 1
        if h == -1:
            h = 590923713
        return h


class TreeNode:
    """
    TreeNode represents a node in a phylogenetic tree.
    Can be either a leaf node (with an index) or intermediate node (with children).
    """
    
    _index: int
    _distance: float
    _is_root: bool
    _parent: Optional[TreeNode]
    _children: List[TreeNode]
    
    def __init__(self, children: Optional[List[TreeNode]] = None, distances: Optional[List[float]] = None, index: Optional[int] = None):
        self._is_root = False
        self._distance = 0.0
        self._parent = None
        
        if index is None:
            # Node is intermediate -> has children
            if children is None or distances is None:
                raise TypeError(
                    "Either reference index (for terminal node) or "
                    "child nodes including the distance "
                    "(for intermediate node) must be set"
                )
            
            if len(children) == 0:
                raise ValueError(
                    "Intermediate nodes must at least contain one child node"
                )
            if len(children) != len(distances):
                raise ValueError(
                    "The number of children must equal the number of distances"
                )
            
            # Check for duplicate children
            for i in range(len(children)):
                for j in range(len(children)):
                    if i != j and children[i] is children[j]:
                        raise ValueError(
                            "Two child nodes cannot be the same object"
                        )
            
            self._index = -1
            self._children = [i for i in children]
            for child, distance in zip(children, distances):
                child._set_parent(self, float(distance))
        elif index < 0:
            raise ValueError("Index cannot be negative")
        else:
            # Node is terminal -> has no children
            if children is not None or distances is not None:
                raise TypeError(
                    "Reference index and child nodes are mutually exclusive"
                )
            self._index = index
            self._children = []
    
    def _set_parent(self, parent: TreeNode, distance: float):
        if self._parent is not None or self._is_root:
            raise ValueError("Node already has a parent")
        self._parent = parent
        self._distance = distance
    
    def copy(self) -> TreeNode:
        """
        Create a deep copy of this TreeNode.
        """
        if self.is_leaf():
            return TreeNode(index=self._index)
        else:
            distances = [child.distance for child in self._children]
            children_clones = [child.copy() for child in self._children]
            return TreeNode(children_clones, distances)
    
    @property
    def index(self) -> Optional[int]:
        return None if self._index == -1 else self._index
    
    @property
    def children(self) -> Optional[List[TreeNode]]:
        return self._children if len(self._children) > 0 else None
    
    @property
    def parent(self) -> Optional[TreeNode]:
        return self._parent
    
    @property
    def distance(self) -> Optional[float]:
        return None if self._parent is None else self._distance
    
    def is_leaf(self) -> bool:
        """Check if the node is a leaf node."""
        return self._index != -1
    
    def is_root(self) -> bool:
        """Check if the node is a root node."""
        return self._is_root
    
    def as_root(self):
        """Convert the node into a root node."""
        if self._parent is not None:
            raise ValueError("Node has parent, cannot be a root node")
        self._is_root = True
    
    def distance_to(self, node: TreeNode, topological: bool = False) -> float:
        """
        Get the distance of this node to another node.
        """
        distance = 0.0
        lca = self.lowest_common_ancestor(node)
        if lca is None:
            raise ValueError("The nodes do not have a common ancestor")
        
        current_node = self
        while current_node is not lca:
            if topological:
                distance += 1.0
            else:
                distance += current_node._distance
            current_node = current_node._parent
        
        current_node = node
        while current_node is not lca:
            if topological:
                distance += 1.0
            else:
                distance += current_node._distance
            current_node = current_node._parent
        
        return distance
    
    def lowest_common_ancestor(self, node: TreeNode) -> Optional[TreeNode]:
        """
        Get the lowest common ancestor of this node and another node.
        """
        lca = None
        # Create two paths from the leaves to root
        self_path = _create_path_to_root(self)
        other_path = _create_path_to_root(node)
        
        # Reverse iteration through path (beginning from root)
        min_len = min(len(self_path), len(other_path))
        for i in range(1, min_len + 1):
            if self_path[-i] is other_path[-i]:
                lca = self_path[-i]
            else:
                break
        return lca
    
    def get_indices(self) -> ndarray:
        """
        Get an array of reference indices that leaf nodes of this node contain.
        """
        return np.array([leaf._index for leaf in self.get_leaves()], dtype=np.int32)
    
    def get_leaves(self) -> List[TreeNode]:
        """
        Get a list of leaf nodes that are direct or indirect child nodes.
        """
        leaf_list = []
        _get_leaves(self, leaf_list)
        return leaf_list
    
    def get_leaf_count(self) -> int:
        """
        Get the number of direct or indirect leaves of this node.
        """
        return _get_leaf_count(self)
    
    def to_newick(self, labels=None, include_distance: bool = True, 
                round_distance: Optional[int] = None) -> str:
        """
        Obtain the node represented in Newick notation.
        """
        if self.is_leaf():
            if labels is not None:
                label = labels[self._index]
                # Check for illegal characters
                illegal_chars = [",", ":", ";", "(", ")"]
                for char in illegal_chars:
                    if char in label:
                        raise ValueError(
                            "Label '" + label + "' contains illegal character '" + char + "'"
                        )
            else:
                label = str(self._index)
            
            if include_distance:
                if round_distance is None:
                    return label + ":" + str(self._distance)
                else:
                    # Use round() and string formatting
                    rounded = round(self._distance, round_distance)
                    return label + ":" + str(rounded)
            else:
                return label
        else:
            # Build string recursively
            child_strings = [child.to_newick(labels, include_distance, round_distance) 
                        for child in self._children]
            child_part = "(" + ",".join(child_strings) + ")"
            
            if include_distance:
                if round_distance is None:
                    return child_part + ":" + str(self._distance)
                else:
                    rounded = round(self._distance, round_distance)
                    return child_part + ":" + str(rounded)
            else:
                return child_part
    @staticmethod
    def from_newick(newick: str, labels: Optional[List[str]] = None) -> Tuple[TreeNode, float]:
        """
        Create a node and all its child nodes from a Newick notation.
        """
        # Ignore any whitespace
        newick = "".join(newick.split())
        
        # Find brackets belonging to sub-newick
        subnewick_start_i = -1
        subnewick_stop_i = -1
        
        for i in range(len(newick)):
            char = newick[i]
            if char == "(":
                subnewick_start_i = i
                break
            if char == ")":
                raise ValueError("Bracket closed before it was opened")
        
        for i in range(len(newick) - 1, -1, -1):
            char = newick[i]
            if char == ")":
                subnewick_stop_i = i + 1
                break
            if char == "(":
                raise ValueError("Bracket was opened but not closed")
        
        if subnewick_start_i == -1 and subnewick_stop_i == -1:
            # No brackets -> Leaf node
            label_and_distance = newick
            try:
                parts = label_and_distance.split(":")
                label = parts[0]
                distance = float(parts[1])
            except:
                distance = 0.0
                label = label_and_distance
            
            if labels is None:
                index = int(label)
            else:
                index = labels.index(label)
            return TreeNode(index=index), distance
        else:
            # Intermediate node
            if subnewick_stop_i == len(newick):
                distance = 0.0
            else:
                label_and_distance = newick[subnewick_stop_i:]
                try:
                    parts = label_and_distance.split(":")
                    distance = float(parts[1])
                except:
                    distance = 0.0
            
            subnewick = newick[subnewick_start_i + 1 : subnewick_stop_i - 1]
            if len(subnewick) == 0:
                raise ValueError("Intermediate node must at least have one child")
            
            # Parse children - split at ',' if at current level
            comma_pos = []
            level = 0
            for i, char in enumerate(subnewick):
                if char == "(":
                    level += 1
                elif char == ")":
                    level -= 1
                elif char == ",":
                    if level == 0:
                        comma_pos.append(i)
                if level < 0:
                    raise ValueError("Bracket closed before it was opened")
            
            children = []
            distances = []
            
            # Recursive tree construction
            for i, pos in enumerate(comma_pos):
                if i == 0:
                    child, dist = TreeNode.from_newick(subnewick[:pos], labels=labels)
                else:
                    prev_pos = comma_pos[i - 1]
                    child, dist = TreeNode.from_newick(subnewick[prev_pos + 1 : pos], labels=labels)
                children.append(child)
                distances.append(dist)
            
            # Node after last comma
            if len(comma_pos) != 0:
                child, dist = TreeNode.from_newick(subnewick[comma_pos[-1] + 1:], labels=labels)
            else:
                child, dist = TreeNode.from_newick(subnewick, labels=labels)
            children.append(child)
            distances.append(dist)
            
            return TreeNode(children, distances), distance
    
    def __str__(self) -> str:
        return self.to_newick()
    
    def __eq__(self, other) -> bool:
        if not isinstance(other, TreeNode):
            return False
        if self._distance != other._distance:
            return False
        if self._index != -1:
            if self._index != other._index:
                return False
        else:
            # Compare children directly without using sets
            if len(self._children) != len(other._children):
                return False
            # Check if all children in self are in other (order doesn't matter)
            for child in self._children:
                found = False
                for other_child in other._children:
                    if child is other_child:  # Use identity comparison
                        found = True
                        break
                if not found:
                    return False
        return True
    def __hash__(self) -> int:
        children_set = set(self._children) if len(self._children) > 0 else set()
        return hash((self._index, children_set, self._distance))


class Tree:
    """
    A Tree represents a rooted phylogenetic tree.
    """
    
    _root: TreeNode
    _leaves: List[TreeNode]
    
    def __init__(self, root: TreeNode):
        root.as_root()
        self._root = root
        
        leaves_unsorted = self._root.get_leaves()
        leaf_count = len(leaves_unsorted)
        
        # Build a temporary dictionary to map indices to leaves
        leaf_dict = {}
        for leaf in leaves_unsorted:
            index = leaf.index
            if index is None:
                raise ValueError("Leaf node must have an index")
            if index >= leaf_count or index < 0:
                raise ValueError("The tree's indices are out of range")
            leaf_dict[index] = leaf
        
        # Build the sorted list
        self._leaves = []
        for i in range(leaf_count):
            if i not in leaf_dict:
                raise ValueError(f"Missing leaf with index {i}")
            self._leaves.append(leaf_dict[i])
    
    def copy(self) -> Tree:
        """Create a deep copy of the tree."""
        return Tree(self._root.copy())
    
    @property
    def root(self) -> TreeNode:
        return self._root
    
    @property
    def leaves(self) -> List[TreeNode]:
        return list(self._leaves)
    
    def get_distance(self, index1: int, index2: int, topological: bool = False) -> float:
        """
        Get the distance between two leaf nodes.
        """
        return self._leaves[index1].distance_to(self._leaves[index2], topological)
    
    def to_newick(self, labels=None, include_distance: bool = True, 
                  round_distance: Optional[int] = None) -> str:
        """
        Obtain the Newick notation of the tree.
        """
        return self._root.to_newick(labels, include_distance, round_distance) + ";"
    
    @staticmethod
    def from_newick(newick: str, labels: Optional[List[str]] = None) -> Tree:
        """
        Create a tree from a Newick notation.
        """
        newick = newick.strip()
        if len(newick) == 0:
            raise ValueError("Newick string is empty")
        # Remove terminal semicolon
        if newick[-1] == ";":
            newick = newick[:-1]
        root, distance = TreeNode.from_newick(newick, labels)
        return Tree(root)
    
    def __str__(self) -> str:
        return self.to_newick()
    
    def __len__(self) -> int:
        return len(self._leaves)
    
    def __eq__(self, other) -> bool:
        if not isinstance(other, Tree):
            return False
        return self._root == other._root
    
    def __hash__(self) -> int:
        return hash(self._root)


def _get_leaves(node: TreeNode, leaf_list: List[TreeNode]):
    """Helper function to recursively collect leaf nodes."""
    if node._index == -1:
        # Intermediate node -> Recursive calls
        for child in node._children:
            _get_leaves(child, leaf_list)
    else:
        # Leaf node -> add to list
        leaf_list.append(node)


def _get_leaf_count(node: TreeNode) -> int:
    """Helper function to count leaf nodes."""
    if node._index == -1:
        # Intermediate node -> Recursive calls
        count = 0
        for child in node._children:
            count += _get_leaf_count(child)
        return count
    else:
        # Leaf node -> return 1
        return 1


def _create_path_to_root(node: TreeNode) -> List[TreeNode]:
    """
    Create a list of nodes representing the path from this node to root.
    """
    path = []
    current_node = node
    path.append(current_node)
    
    while current_node._parent is not None:
        current_node = current_node._parent
        path.append(current_node)
    
    return path


def as_binary(tree_or_node):
    """
    Convert a tree into a binary tree.
    """
    if isinstance(tree_or_node, Tree):
        node, _ = _as_binary(tree_or_node.root)
        return Tree(node)
    elif isinstance(tree_or_node, TreeNode):
        node, _ = _as_binary(tree_or_node)
        return node
    else:
        raise TypeError(
            f"Expected 'Tree' or 'TreeNode', not {type(tree_or_node).__name__}"
        )


def _as_binary(node: TreeNode) -> Tuple[TreeNode, Optional[float]]:
    """
    The actual logic for converting to binary tree.
    """
    children = node.children
    
    if children is None:
        # Leaf node
        return TreeNode(index=node.index), node.distance
    elif len(children) == 1:
        # Intermediate node with one child -> omit node
        child, distance = _as_binary(node.children[0])
        if node.is_root():
            return child, None
        else:
            return child, node.distance + distance
    elif len(children) > 2:
        # More than two children -> create binary nodes
        rem_children = []
        distances = []
        for child in children:
            bin_child, dist = _as_binary(child)
            rem_children.append(bin_child)
            distances.append(dist)
        
        current_div_node = None
        while len(rem_children) > 0:
            if current_div_node is None:
                # Bottom-most node with two children
                current_div_node = TreeNode(
                    rem_children[:2],
                    distances[:2]
                )
                rem_children = rem_children[2:]
                distances = distances[2:]
            else:
                # Node with one remaining child and previous intermediate
                current_div_node = TreeNode(
                    [current_div_node, rem_children[0]],
                    [0.0, distances[0]]
                )
                rem_children = rem_children[1:]
                distances = distances[1:]
        
        return current_div_node, node.distance
    else:
        # Exactly two children -> keep unchanged
        binary_children = []
        bin_distances = []
        for child in children:
            bin_child, dist = _as_binary(child)
            binary_children.append(bin_child)
            bin_distances.append(dist)
        return TreeNode(binary_children, bin_distances), node.distance