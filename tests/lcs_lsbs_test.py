#!/usr/bin/env python3
"""
Unit tests for the LSBS-based LCS implementation.

Copyright (C) 2025 Awwal Mohammed, Caroline Sumathi Selvarajah

This software is dual-licensed under MIT OR GPL-3.0.
Choose the license that best fits your use case.
See LICENSE file for full license terms.
"""
import os,sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import unittest

# noinspection PyPep8Naming
from lcs_lsbs import lcs_lsbs
from lis_funcs import *


class MyTestCase(unittest.TestCase):

    def setUp(self):
        self.pairs = [
            ("ABCBDAB", "BDCAB"),  # BDAB (n=4)
            ("GTTCCTAATA", "CGATAATTGAGA"),  # GTTTAA or GTAATA (6)
            ("GTTCCTAAT", "CGATAATTGAG"),  # GTTTA or GTAAT (5)
            ("12345", "23415"),  # 2345 (4)
            ("123", "1323"),  # 123 (3)
            ("126548", "216544"),  # 2654 or 1654 (4)
            ("AGGTAB", "GXTXAYB"),  # GTAB (4)
            ("BD", "ABCD"),  # BD (2)
            ("ABCDGH", "AEDFHR"),  # ADH (3)
            ("ABCDE", "ACE"),  # ACE (3)
            ("hofubmnylkra", "pqhgxgdofcvmr"),  # hofmr (5)
            ("oxcpqrsvwf", "shmtulqrypy"),  # "qr"
            ("ABC", "DEF"),  # ∅ (0)
            ("XYYYYXY", "XYYYYY"),  # XYYYYY
        ]

        self.expected_results = [
            (4, "BDAB"),
            (6, "GTTTAA"),
            (5, "GTTTA"),
            (4, "2345"),
            (3, "123"),
            (4, "2654"),  # or "1654" - this could vary
            (4, "GTAB"),
            (2, "BD"),
            (3, "ADH"),
            (3, "ACE"),
            (5, "hofmr"),
            (2, "qr"),  # crucial edge case
            (0, ""),
            (6, "XYYYYY")  # bad quadratic space strings
        ]

    def test_lcs_cet(self):
        for i, (A, B) in enumerate(self.pairs):
            with self.subTest(i=i):
                dp_result_len, dp_result_str = lcs_lsbs(a=A, b=B, f=lis_dp_2d)
                ps_result_len, ps_result_str = lcs_lsbs(a=A, b=B, f=lis_ps_2d)
                expected_len, expected_str = self.expected_results[i]

                # Print the expected and actual strings
                print(f"Test Case {i + 1}:")
                print(f"A = {A}  B = {B}")
                print(f"Expected String: {expected_str}")
                print(f"LCS via DP LIS :: Actual String: {dp_result_str}")
                print(f"LCS via PS LIS :: Actual String: {ps_result_str}")
                print()

                # Test only for the length
                self.assertEqual(dp_result_len, expected_len)
                self.assertEqual(ps_result_len, expected_len)

        # self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
