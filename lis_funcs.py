#!/usr/bin/env python3
"""
A collection of pair-wise variants of the Longest Increasing Subsequence (LIS) algorithm.
- Classic Dynamic Programming LIS, 2D variant.
- Patience Sorting LIS, 1D & 2D variants.

Copyright (C) 2025 Awwal Mohammed, Caroline Sumathi Selvarajah

This software is dual-licensed under MIT OR GPL-3.0.
Choose the license that best fits your use case.
See LICENSE file for full license terms.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple
from bisect import bisect_left


def lis_ps_2d(nums: List[Tuple[int, int]]):
    """
    Finds the Longest Increasing Subsequence for pairs (value, bucket_index),
    where both value and bucket_index must be strictly increasing.
    Uses the standard O(n Log n) algorithm for 2D LIS:
    1. Sort points (value, bucket_index) by value asc, then bucket_index desc.
    2. Find LIS on the bucket_index components.

    Identical overriding technique used herewith is explained in this video:
    https://www.youtube.com/watch?v=OIU8ZLC4qIQ
    """
    n = len(nums)
    if n == 0:
        return []

    # Enrich nums with original indices to reconstruct the actual values later if needed,
    # though this specific function returns the 'value' (first component) of the LIS pairs.
    # Format: (value, bucket_index, original_index_in_nums)
    # For this problem, nums are (b_idx, a_idx). We need b_idx_i < b_idx_j and a_idx_i < a_idx_j.
    # So, x = p[0] (b_idx), y = p[1] (a_idx)
    
    # Create a list of tuples: (value, bucket_index, original_item_tuple)
    # We need to return List[int] of the first components (values) of the LIS pairs.
    # The original nums list is List[Tuple[value, bucket_index]]
    
    # Sort by value (nums[i][0]) ascending, then bucket_index (nums[i][1]) descending.
    # Store original pair to reconstruct LIS of pairs.
    # sorted_nums_enriched will store (value, bucket_index, original_pair_from_nums)
    # No, we need original *index* if we want to reconstruct LIS from original nums,
    # but the problem asks for a list of values.
    # The LIS construction here will give us the LIS of (value, bucket_index) pairs.
    
    # Create (value, bucket_index) items. The LIS function will operate on these.
    # The key for sorting: item[0] is value, item[1] is bucket_index.
    # Sort by item[0] asc, item[1] desc.
    # We are passing `nums` directly. The elements are (b_idx, a_idx).
    # So, sort by b_idx (nums[i][0]) asc, then a_idx (nums[i][1]) desc.
    
    # Create a list of items to sort, preserving original pairs for final LIS.
    # Each element: (value, bucket_index, original_pair_tuple)
    # Using original_pair_tuple is fine as we extract value from it.

    # ADDITIONAL NOTES ON THIS BILATERAL SORTING MECHANISM
    # .sort(key=lambda x: (x[0], -x[1])) means the concurrent lowest set x[0]'s
    # are taken as a group in ASCENDING order, then each correspionding x[1]'s
    # are ordered for each group in DESCENDING order.
    enriched_items = [(item[0], item[1], item) for item in nums]
    ei2 = enriched_items.copy()
    enriched_items.sort(key=lambda x: (x[0], -x[1]))

    if not enriched_items: # Should be caught by n==0 but good practice
        return []

    # Now find LIS on the bucket_index components (item[1]) of enriched_items
    # Standard O(N log N) LIS algorithm (e.g., using patience sorting idea)
    # `tails_values` stores the smallest ending bucket_index of an increasing subsequence of buckets of length k+1.
    tails_values = []
    # `lis_candidate_pairs` stores the actual (value, bucket_index) pair from `enriched_items`
    # that corresponds to the entry in `tails_values`. This helps reconstruct the LIS of pairs.
    # More precisely, `tails_end_pairs[k]` is the pair ending an LIS of length k+1.
    tails_end_pairs: List[Tuple[int,int]] = [] # Stores the (value, bucket_index) pair, not just bucket

    # `predecessors_for_lis_pairs` maps an enriched_item (by identity or index)
    # to its predecessor enriched_item in the LIS. For simplicity, we store the pair itself.
    # To reconstruct, it's better to store indices or use a standard LIS reconstruction with P array.

    # Simpler LIS reconstruction:
    # M[j] stores the index in `enriched_items` of the smallest tail of all LIS of buckets of length j.
    # preds[i] stores the index in `enriched_items` of the predecessor of enriched_items[i].
    
    num_enriched = len(enriched_items)
    preds = [-1] * num_enriched # Predecessor indices in enriched_items
    
    # M[k] stores index in enriched_items for LIS of length k (1-indexed length)
    # To map to 0-indexed tails_values: M[k] is index for tails_values[k-1]
    # Let active_lis_tails_positions[k] be index in enriched_items of item ending LIS of length k+1
    active_lis_tails_positions = []
    item_values = []

    for i in range(num_enriched):
        current_item_value = enriched_items[i][0] # This is original value (b_idx)
        current_item_bucket = enriched_items[i][1] # This is original bucket (a_idx)

        # Find insertion point for current_item_bucket in tails_values
        # `j` is the length of LIS ending with predecessor + 1 (0-indexed)
        # i.e., current_item_bucket would be the end of an LIS of length `j+1`.
        j = bisect_left(tails_values, current_item_bucket)

        if j == len(tails_values):
            tails_values.append(current_item_bucket) # Store bucket index of enriched_items pair
            active_lis_tails_positions.append(i) # Store index from enriched_items
        else:
            tails_values[j] = current_item_bucket  # Override bucket_index of enriched_items pair
            active_lis_tails_positions[j] = i # Override index from enriched_items
        
        # Set predecessor
        # The item at enriched_items[i] extends the LIS ending at
        # enriched_items[active_lis_tails_positions[j-1]] if j > 0.
        if j > 0:
            preds[i] = active_lis_tails_positions[j-1]
        # else preds[i] remains -1 (starts an LIS of length 1)

    # Reconstruct LIS of (value, bucket_index) pairs
    lis_of_pairs = []
    if active_lis_tails_positions: # If any LIS was found
        # Start from the end of the longest LIS found
        # The last element of `active_lis_tails_indices` is the index (in `enriched_items`)
        # of the item that ends the overall LIS.
        current_idx_in_enriched = active_lis_tails_positions[-1]
        while current_idx_in_enriched != -1:
            # The actual pair is enriched_items[current_idx_in_enriched][2] (the original tuple)
            # or simply enriched_items[current_idx_in_enriched] if we only stored (value, bucket)
            # We need the value component of this pair for the final result.
            # The pair itself is (value, bucket_index, original_item_from_nums)
            # We need the value (enriched_items[current_idx_in_enriched][0])
            lis_of_pairs.append(enriched_items[current_idx_in_enriched][2]) # Append the original (value, bucket) tuple
            current_idx_in_enriched = preds[current_idx_in_enriched]
        lis_of_pairs.reverse()

    # print(
    #     "CHECKING TAILS AND LIS CONTENTS...", 
    #     [p[0] for p in lis_of_pairs], sep="\n"
    # )
    return [p[0] for p in lis_of_pairs] # Return first components (values)


def lis_dp_2d(nums: List[Tuple[int, int]]):
    """
    Classic dynamic programming LIS solution adapted for 2D/pair-wise bucket constraint.
    Takes O(n Log n) time, not as efficient as Patience Sorting counterpart.
    """
    n = len(nums)
    dp = [1] * n                            # Assume each element is the singleton LIS
    prev = [-1] * n                         # Previous index in the longest sequence

    longest_seq_end = 0
    for i in range(1, n):
        for j in range(i):
            if (
                nums[j][1] < nums[i][1] and   # Invariant 1: Check if jth is from an earlier bucket than ith
                nums[j][0] < nums[i][0] and   # Invariant 2: With ith number always ahead, jth MUST always be less
                dp[j] + 1 > dp[i]             # Invariant 3: Only consider iterations that INCREASE the LIS
            ):                     
                dp[i] = dp[j] + 1             # Iteration i yields a LONGER subsequence (or update if better candidate is found)
                prev[i] = j                   # j is a direct valid LIS predecessor of i, so add a back pointer

        if dp[i] > dp[longest_seq_end]:       # Determine position of longest sequence so far
            longest_seq_end = i

    # return longest_seq_end                  # If we want only length of LCS

    sequence = []                             # To hold reconstructed sequence
    if len(nums) > 0:
        while longest_seq_end != -1:
            sequence.append(nums[longest_seq_end][0])
            longest_seq_end = prev[longest_seq_end]

    return sequence[::-1]                     # Reverse to get the sequence in ascending order


def lis_ps_1d(nums: List[int]):
    """
    Constructs the (1-Dimensional) LIS for a sequence of integers by segregating monotically sorted piles
    of integers, whose top elements are indexed in a (prefrentially implicit) binary search tree,
    as done in the Patience Sorting algorithm.
    """
    @dataclass
    class Node:
        stack: List[int]
        top: int
        nxt: Node | None
        candidate: int | None

    bst: List[Node] = []  # Pseudo-BST (only search ops are logarithmic)
    n = len(nums)

    for i in range(n):
        ins_point = bisect_left(bst, nums[i], key=lambda n: n.top)
        if len(bst) == 0:
            bst.insert(ins_point, Node([nums[i]], nums[i], None, nums[i]))
        else:
            if ins_point < len(bst):
                if bst[ins_point].top > nums[i]:
                    bst[ins_point].stack.append(nums[i])
                    bst[ins_point].top = nums[i]
            else:   # Create a new pile, i.e. insert new BST Node, and add pointer to LIS candidate
                bst[ins_point-1].candidate = bst[ins_point-1].top
                new_node = Node([nums[i]], nums[i], None, nums[i])
                bst.insert(ins_point, new_node)
                bst[ins_point-1].nxt = new_node 

    return [n.candidate for n in bst]   


#-------------------------------------------------------------------
def test(inputs: List[List[int]]):
    print("Testing 1D LIS using Patience Sorting...")
    for i in inputs:
        print(lis_ps_1d(i))


if __name__ == '__main__':
    data = [
        [1,3,2,3],
        [10, 22, 9, 33, 21, 50, 41, 60],
        [1, 2, 3, 4],
        [1, 3, 2, 4],
        [3, 10, 2, 1, 20],
        [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15],
        [5, 1, 6, 2, 7, 3, 8]
    ]
    test(data)
    