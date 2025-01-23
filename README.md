# RNAProcessor

A comprehensive Python pipeline for RNA sequence processing and BLAST analysis.

## Overview

RNAProcessor is a modular tool designed to process RNA sequences through multiple filtering and analysis steps. The pipeline includes duplicate removal, adapter trimming, quality filtering, pair matching, alignment analysis, and BLAST searching.

## Pipeline Steps and Logic

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
     - Generates formatted PDF reports with:
       - Hit summaries
       - Alignment details
       - Quality metrics

## Implementation Details

### Parallel Processing
- Uses ProcessPoolExecutor for CPU-intensive tasks
- Configurable core count (default: 22)
- Chunk-based processing for memory efficiency

### Error Handling
- Comprehensive logging at each step
- Detailed error messages and stack traces
- Statistics generation for quality control

### File Management
- Organized directory structure
- Consistent file naming conventions
- Intermediate results preserved for analysis

## Installation
