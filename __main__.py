#!/usr/bin/env python3
"""
Interactive mode for executing LCS via LSBS

Copyright (C) 2025 Awwal Mohammed, Caroline Sumathi Selvarajah

This software is dual-licensed under MIT OR GPL-3.0.
Choose the license that best fits your use case.
See LICENSE file for full license terms.
"""

from lcs_lsbs import lcs_lsbs


def main():
    while True:
        print("Press CTRL+C to quit\n")
        x = input("Enter first sequence:\t")
        y = input("Enter second sequence: \t")
        print(lcs_lsbs(x, y))


if __name__ == "__main__":
    main()