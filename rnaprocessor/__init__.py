"""
RNAProcessor - A module for RNA sequence processing and BLAST analysis
"""

from rnaprocessor.src.duplicate_remover import DuplicateRemover
from rnaprocessor.src.adapter_remover import AdapterRemover
from rnaprocessor.src.quality_filter import QualityFilter
from rnaprocessor.src.pair_matcher import PairMatcher
from rnaprocessor.src.alignment_analyzer import AlignmentAnalyzer
from rnaprocessor.src.blast_searcher import BlastSearcher
from rnaprocessor.src.pdf_generator import BlastPDFGenerator

__version__ = "0.1.0"
__author__ = "RNA Analysis Team"

# Export all components for easy access
__all__ = [
    'DuplicateRemover',
    'AdapterRemover',
    'QualityFilter',
    'PairMatcher',
    'AlignmentAnalyzer',
    'BlastSearcher',
    'BlastPDFGenerator',
] 