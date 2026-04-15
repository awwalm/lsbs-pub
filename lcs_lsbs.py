#!/usr/bin/env python3
"""
Longest Common Subsequence using our own Longest Sorted Bucket Sequence intuition.

Copyright (C) 2025 Awwal Mohammed, Caroline Sumathi Selvarajah

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from cet import get_cet, improved_cet
from lsbs import lsbs
from typing import Sequence, Callable


def lcs_lsbs(a: Sequence, b: Sequence, f: Callable = None):
    """
    Compute Longest Common Subsequence based on CET,
    using the Longest Sorted Bucket Sequence (LSBS) algorithm.
    :param f: A desired LIS callback algorithm on the matching buckets triggered by LSBS.
    """
    # matching_buckets = get_cet(a, b)[1]         # We only need the queues (second row of CET)
    matching_buckets = improved_cet(a, b)
    indices = lsbs(buckets=matching_buckets) if not f else lsbs(buckets=matching_buckets, lis_func=f)
    # print("indices from lsbs", indices)
    matches = [b[i] for i in indices]
    return len(matches), str().join(matches)
