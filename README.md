# Longest Sorted Bucket Sequence: An Algorithmic Bridge Between LCS and 2D LIS

[![License: Dual](https://img.shields.io/badge/License-MIT%20OR%20GPL--3.0-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

Official implementation and research repository for the **Longest Sorted Bucket Sequence (LSBS)** algorithm, a unique algorithmic approach that bridges the Longest Common Subsequence (LCS) and 2D Longest Increasing Subsequence (LIS) problems.

## 🔬 Overview

LSBS is a unique algorithm that solves the Longest Common Subsequence problem by transforming it into a 2D Longest Increasing Subsequence problem through the use of a **Common Element Table (CET)**. The approach provides both theoretical insights and practical implementations for sequence analysis.

### Key Innovation

- **LSBS Algorithm**: Transforms LCS computation into 2D LIS finding
- **Common Element Table (CET)**: Efficiently maps character occurrences between sequences
- **Multiple LIS Backends**: Supports both Dynamic Programming and Patience Sorting approaches
- **Performance Analysis**: Comprehensive comparison with classical LCS algorithms

## 🚀 Quick Start

### Interactive Mode
```bash
python __main__.py
```

### Programmatic Usage
```python
from lcs_lsbs import lcs_lsbs

# Basic usage
sequence1 = "ABCBDAB"
sequence2 = "BDCAB"
length, lcs = lcs_lsbs(sequence1, sequence2)
print(f"LCS length: {length}, LCS: '{lcs}'")  # Output: LCS length: 4, LCS: 'BDAB'

# Using different LIS algorithms
from lis_funcs import lis_ps_2d
length, lcs = lcs_lsbs(sequence1, sequence2, f=lis_ps_2d)
```

## 📊 Algorithm Implementations

This repository includes multiple algorithm implementations for comparative analysis:

### LCS Algorithms
- **LCS-LSBS** (`lcs_lsbs.py`): Our novel LSBS-based approach
- **Classic DP** (`lcs_dp.py`): Traditional dynamic programming solution
- **Recursive variants** (`lcs_rec.py`): Naive, memoized, and cached implementations

### LIS Algorithms
- **2D DP LIS** (`lis_funcs.py`): Classic O(n²) dynamic programming
- **2D Patience Sorting** (`lis_funcs.py`): Optimized O(n log n) approach
- **1D Patience Sorting** (`lis_funcs.py`): Traditional LIS with patience sorting

### Supporting Components
- **Common Element Table** (`cet.py`): Core data structure for LSBS
- **LSBS Core** (`lsbs.py`): Main LSBS algorithm implementation

## 🏗️ Project Structure

```
.
├── __main__.py           # Interactive CLI interface
├── lcs_lsbs.py          # Main LSBS-based LCS implementation
├── lsbs.py              # Core LSBS algorithm
├── cet.py               # Common Element Table implementation
├── lis_funcs.py         # Collection of LIS algorithm variants
├── lcs_dp.py            # Classic DP LCS implementation
├── lcs_rec.py           # Recursive LCS implementations
├── tests/
│   └── lcs_lsbs_test.py # Unit tests and benchmarks
├── LICENSE              # GPL v3 License
└── README.md            # This file
```

## 🧪 Testing

Run the comprehensive test suite:

```bash
python -m pytest tests/
```

Or run individual algorithm tests:

```bash
# Test LSBS implementation
python lsbs.py

# Test classic DP LCS
python lcs_dp.py

# Test recursive LCS variants
python lcs_rec.py

# Test LIS algorithms
python lis_funcs.py
```

## 📈 Performance Characteristics

| Algorithm | Time Complexity | Space Complexity | Notes |
|-----------|----------------|------------------|-------|
| LSBS (DP backend) | O(n + m + k²) | O(nm) | k = number of matches |
| LSBS (PS backend) | O(n + m + k log k) | O(nm) | Optimal for sparse matches |
| Classic DP LCS | O(nm) | O(nm) | Traditional approach |
| Recursive LCS | O(2^(n+m)) | O(n+m) | Exponential without memoization |

## 🔧 Algorithm Details

### Common Element Table (CET)
The CET maps each character in sequence A to all its occurrence positions in sequence B:
- **Row 1**: Indices of characters in sequence A
- **Row 2**: Queues of matching positions in sequence B

### LSBS Process
1. **Build CET**: Create mapping of character occurrences
2. **Flatten to pairs**: Convert to (`position_in_B, bucket_index`) pairs
3. **Apply 2D LIS**: Find longest increasing subsequence on both dimensions
4. **Reconstruct LCS**: Map result back to original characters

## 📚 Research Context

This implementation accompanies the research manuscript:
> "Longest Sorted Bucket Sequence: An Algorithmic Bridge Between LCS and 2D LIS"

The work demonstrates how classical sequence problems can be viewed through different algorithmic lenses, potentially offering new insights for:
- Bioinformatics sequence alignment
- Text similarity analysis
- Pattern matching applications

## 🛠️ Requirements

- Python 3.8+
- No external dependencies (uses only standard library)

## 📄 License

This project is **dual-licensed** to provide maximum flexibility:

### Choose Your License:

🔓 **MIT License** (More Permissive)
- ✅ Commercial use without restrictions
- ✅ Can be included in proprietary software  
- ✅ Minimal requirements (just attribution)
- 👥 **Ideal for**: Companies, proprietary projects, maximum adoption

🔒 **GPL 3.0 License** (Copyleft)
- ✅ Commercial use allowed
- ⚠️ Derivative works must also be open source
- ⚠️ Must provide source code when distributing
- 🌍 **Ideal for**: Open source projects, research collaboration, ensuring code stays open

### How to Choose:
- **For commercial/proprietary use**: Choose MIT License
- **For open source projects**: Choose GPL 3.0 License  
- **For research**: Either license works, GPL 3.0 ensures derivatives stay open

See the [LICENSE](LICENSE) file for complete license texts.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## 📧 Citation

If you use this work in your research, please cite:
```bibtex
@software{lsbs2025,
  title={Longest Sorted Bucket Sequence: An Algorithmic Bridge Between LCS and 2D LIS},
  author={[Awwal Mohammed; Caroline Sumathi Selvarajah]},
  year={2025},
  url={https://github.com/awwalm/lsbs-public}
}
```

---

**Note**: This is a research implementation. For production use, consider performance optimizations and additional testing for your specific use case.
