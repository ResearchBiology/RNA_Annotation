# RNAProcessor

A comprehensive Python pipeline for RNA sequence processing and BLAST analysis.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Steps](#pipeline-steps)
- [Usage](#usage)
- [Implementation Details](#implementation-details)
- [Contributing](#contributing)

## Overview

RNAProcessor is a modular tool designed to process RNA sequences through multiple filtering and analysis steps. The pipeline includes duplicate removal, adapter trimming, quality filtering, pair matching, alignment analysis, and BLAST searching.

## Installation

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Steps
1. Clone the repository:
```bash
git clone https://github.com/Cheatallin/rnaprocessor.git
cd rnaprocessor
```

2. Install the package and its dependencies:
```bash
pip install .
```

Required dependencies will be automatically installed:
- biopython >= 1.80
- numpy >= 1.21.0
- pandas >= 1.3.0
- matplotlib >= 3.4.0
- seaborn >= 0.11.0

## Quick Start

1. Prepare your input files:
   - Place your `.fq` files in an input directory
   - Each file should be properly formatted with:
     - Line 1: Sequence information
     - Line 2: Sequence data

2. Run the pipeline:
```bash
rnaprocessor --input-dir path/to/input/files --output-dir path/to/output
```

The pipeline will automatically:
1. Remove duplicates (initial)
2. Remove adapters
3. Filter by quality
4. Match sequence pairs
5. Remove duplicates (secondary)
6. Analyze alignments
7. Perform BLAST search

## Pipeline Steps

1. **Initial Duplicate Removal** (`DuplicateRemover`)
   - **Purpose**: Remove redundant sequence pairs to reduce computational load
   - **Logic**:
     - Reads sequence pairs from input .fq files
     - Creates tuple keys from sequence pairs
     - Uses Python set for O(1) lookup of duplicates
     - Only writes first occurrence of each unique sequence pair
     - Tracks statistics: total pairs, unique pairs, duplication rate

2. **Adapter Removal** (`AdapterRemover`)
   - **Purpose**: Remove adapter sequences that could interfere with analysis
   - **Logic**:
     - Searches for identifier tag "ATAGCGG" in last 30bp of sequence
     - If found: Trims sequence up to identifier position
     - If not found: Removes entire sequence
     - Maintains sequence integrity by precise trimming
     - Generates detailed statistics about trimming results

3. **Quality Score Filtering** (`QualityFilter`)
   - **Purpose**: Ensure high sequence quality for downstream analysis
   - **Logic**:
     - Uses sliding window approach (default size: 15)
     - Converts Phred+33 ASCII quality scores to numeric values
     - Calculates mean quality score within each window
     - Rejects sequences if any window falls below Q30
     - Generates quality distribution plots for visualization
     - Tracks position-wise quality statistics

4. **Pair Matching** (`PairMatcher`)
   - **Purpose**: Match forward and reverse reads based on overlap
   - **Logic**:
     - Identifies forward/reverse file pairs by naming pattern
     - Uses minimum 30bp overlap requirement
     - Implements parallel processing for large datasets
     - Searches for significant overlapping regions between pairs
     - Validates matches using sequence alignment
     - Generates paired output maintaining read order

5. **Secondary Duplicate Removal** (`DuplicateRemover`)
   - **Purpose**: Remove duplicates that emerge after pairing
   - **Logic**:
     - Similar to initial duplicate removal
     - Operates on paired sequences
     - Considers both sequences in pair for uniqueness
     - Updates statistics post-pairing

6. **Alignment Analysis** (`AlignmentAnalyzer`)
   - **Purpose**: Identify most frequent sequences for BLAST
   - **Logic**:
     - Processes paired sequence files
     - Uses Counter object for frequency tracking
     - Identifies most common sequence per sample
     - Records sequence metadata and statistics
     - Prepares sequences for BLAST analysis

7. **BLAST Search** (`BlastSearcher` and `BlastPDFGenerator`)
   - **Purpose**: Identify sequence origins and generate reports
   - **Logic**:
     - Performs BLASTN search against NT database
     - Implements rate limiting for NCBI API
     - Parses XML BLAST results
     - Extracts key alignment information:
       - E-values
       - Identity percentages
       - Gap statistics

## Usage

### Command Line Arguments

```bash
rnaprocessor [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR] [--cores CORES]
```

#### Arguments:
- `--input-dir`: Directory containing input .fq files (default: "label_cleaned_files")
- `--output-dir`: Base directory for output files (default: "results")
- `--cores`: Number of CPU cores to use (default: 22)

### Output Directory Structure

The pipeline creates the following directory structure for outputs:
```
output_dir/
├── duplicates/        # Initial duplicate removal results
├── adapters/         # Adapter removal results
├── quality/          # Quality filtering results
├── paired/          # Pair matching results
├── new_duplicate/   # Secondary duplicate removal results
├── alignment/       # Alignment analysis results         
└── blast/         # BLAST search results
```

### Progress Tracking

The pipeline provides detailed logging of each step, including:
- Number of files processed
- Success/failure status
- Statistics for each processing step
- Final summary of the entire pipeline

### Error Handling

- The pipeline includes comprehensive error handling
- Each step is logged with detailed error messages
- Processing continues even if individual files fail
- Final summary includes success/failure counts

## Implementation Details

### Parallel Processing
- Uses ProcessPoolExecutor for CPU-intensive tasks
- Default: 22 cores (configurable via --cores)
- Automatically adjusts to available CPU cores

### File Management
- Input files must be in .fq format
- Organized directory structure for each processing step
- Intermediate results are preserved for analysis
- Consistent file naming conventions maintained throughout the pipeline

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.
