#!/usr/bin/env python3
"""
The Patience Sorting algorithm reveals how the classic dynamic programming solution 
to the Longest Increasing Subsequence (LIS) can be reduced from quadratic time to log-linear time.

Copyright (C) 2025 Awwal Mohammed, Caroline Sumathi Selvarajah

This software is dual-licensed under MIT OR GPL-3.0.
Choose the license that best fits your use case.
See LICENSE file for full license terms.
"""
from __future__ import annotations
from bisect import bisect_left
from dataclasses import dataclass
from typing import List

