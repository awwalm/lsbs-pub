#!/usr/bin/env python3
"""
The Longest Sorted Bucket Sequence (LSBS) is a sequence formed by taking a candidate 
from each sorted bucket, without repetitions (if re-encountered in a later bucket) 
in a strictly increasing manner.

Copyright (C) 2025 Awwal Mohammed, Caroline Sumathi Selvarajah

This software is dual-licensed under MIT OR GPL-3.0.
Choose the license that best fits your use case.
See LICENSE file for full license terms.
"""
from typing import List, Callable
from lis_funcs import *


def lsbs(buckets: List[List[int]], lis_func: Callable=lis_dp_2d, sort=False) -> List[int]:
    """
    Extracts the Longest Sorted Bucket Sequence (LSBS) from an array of buckets 
    by adapting the Longest Increasing Subsequence (LIS) algorithm.
    Take note that the buckets computed by the prevailing CET function
    between sequences A and B is sorted with respect due to A.

    :param buckets: The second row (buckets of indices) of the Common Element Table (see cet.py).
    :param lis_func: The LIS callback that takes as argument a list of int pairs (default is DP algorithm).
    :param sort: Optional boolean flag (default is False) for issuing a sort command on the buckets. 
    :returns: The longest sorted sequence of indices of string b, that are characters of string a.
    """
    flattened = []
    for i, b in enumerate(buckets):
        for num in b:
            flattened.append((num, i))      # Append tuple (num, i) to the flattened list
    if sort:
        flattened = sorted(flattened)       # Sort the flattened list if required
    
    return lis_func(flattened)


if __name__ == "__main__":
    # Test cases
    test_cases = [
        [[3], [0, 4], [2], [0, 4], [1], [3], [0, 4]],
        [[1, 8, 10], [3, 6, 7], [3, 6, 7], [0], [0], [3, 6, 7], [2, 4, 5, 9, 11], [2, 4, 5, 9, 11], [3, 6, 7],
         [2, 4, 5, 9, 11]],
        [[1, 8, 10], [3, 6, 7], [3, 6, 7], [0], [0], [3, 6, 7], [2, 4, 5, 9], [2, 4, 5, 9], [3, 6, 7]],
        [[3], [0], [1], [2], [4]],
        [[], [], [], [9], [6], [7], [0], [], [], []]
    ]

    for t, case in enumerate(test_cases, 1):
        result_dp = lsbs(case, lis_dp_2d)
        result_ps = lsbs(case, lis_ps_2d)
        print(f"Test case {t}:")
        print(f"Input: {case}")
        print(f"Output via DP LIS: {result_dp}")
        print(f"Output via PS/BS LIS: {result_ps}")
        print(f"Length: {len(result_dp)}\n")
