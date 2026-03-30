#!/usr/bin/env python3
"""
The Common Element Table (CET) (of two strings) consists of matching pairs of each character
in sequence A, mapped to all its existing positions in sequence B.

Copyright (C) 2025 Awwal Mohammed, Caroline Sumathi Selvarajah

This software is dual-licensed under MIT OR GPL-3.0.
Choose the license that best fits your use case.
See LICENSE file for full license terms.
"""

from collections import deque
from typing import Sequence


def get_cet(a: Sequence, b: Sequence):
    """The Common Element Table function for obtaining
    all indices matching each character in string A to string B.
    """
    cet: list[list] = [                     # Common Element Table
        [j for j in range(len(a))],         # Row 1: Indices of string A
        [deque() for _ in range(len(a))]    # Row 2: Empty queue for occurences of A[j] in B
    ]

    for j in range(len(a)):                 # Construct CET with respect due to A
        for k in range(len(b)):             # Find all occurences of A[j] in B
            if a[j] == b[k]:
                cet[1][j].append(k)         # Record matched indices in queue

    return cet


def improved_cet(a: Sequence, b: Sequence):
    """Same as CET but more optimal and returns only matching buckets.
    - Expected time complexity: O(a + b)
    - Worst case time complexity: O(a * b); arises in 'abnormal' strings.
    """
    # Build character-to-indices mapping for string B in O(m) time
    char_to_indices = {}
    for k in range(len(b)):
        char = b[k]
        if char not in char_to_indices:
            char_to_indices[char] = []
        char_to_indices[char].append(k)
    
    # Build result for each character in A in O(n) time
    matching_buckets = []
    for j in range(len(a)):
        char = a[j]
        if char in char_to_indices:
            # Create a deque with matching indices (to match original structure)
            matching_buckets.append(deque(char_to_indices[char]))
        else:
            # Empty deque if no matches found
            matching_buckets.append(deque())
    
    return matching_buckets
