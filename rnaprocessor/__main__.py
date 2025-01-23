import argparse
import logging
from pathlib import Path
import sys
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

from .src.duplicate_remover import DuplicateRemover
from .src.adapter_remover import AdapterRemover
from .src.quality_filter import QualityFilter
from .src.pair_matcher import PairMatcher
from .src.alignment_analyzer import AlignmentAnalyzer
from .src.blast_searcher import BlastSearcher
from .src.pdf_generator import BlastPDFGenerator

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def get_input_files(input_dir="label_cleaned_files"):
    """Get all .fq files from input directory"""
    input_path = Path(input_dir)
    if not input_path.exists():
        raise ValueError(f"Input directory not found: {input_dir}")
    
    files = list(input_path.glob("*.fq"))
    if not files:
        raise ValueError(f"No .fq files found in {input_dir}")
    
    return files

def main():
    parser = argparse.ArgumentParser(description="RNA Processing Pipeline")
    
    parser.add_argument(
        '--input-dir',
        type=str,
        default="label_cleaned_files",
        help='Directory containing input .fq files'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default="results",
        help='Base directory for output files'
    )
    
    parser.add_argument(
        '--cores',
        type=int,
        default=22,
        help='Number of CPU cores to use'
    )
    
    args = parser.parse_args()
    
    try:
        # Set number of cores for multiprocessing
        num_cores = min(args.cores, multiprocessing.cpu_count())
        multiprocessing.set_start_method('spawn', force=True)
        
        logger.info(f"Starting RNA processing pipeline with {num_cores} cores")
        logger.info("Pipeline steps: Initial Duplicate Removal -> Adapter Removal -> Quality Filtering -> "
                   "Pair Matching -> Secondary Duplicate Removal -> Alignment Analysis -> BLAST Search")
        
        # Get input files
        input_files = get_input_files(args.input_dir)
        logger.info(f"Found {len(input_files)} input files to process")
        
        # 1. Initial Duplicate Removal
        logger.info("\n=== Step 1: Initial Duplicate Removal ===")
        duplicate_remover = DuplicateRemover(output_dir=f"{args.output_dir}/duplicates")
        duplicate_files = []
        for input_file in input_files:
            result = duplicate_remover.remove_duplicates(input_file)
            if result:
                duplicate_files.append(result)
        logger.info(f"Completed initial duplicate removal. Generated {len(duplicate_files)} files")
        
        # 2. Adapter Removal
        logger.info("\n=== Step 2: Adapter Removal ===")
        adapter_remover = AdapterRemover(output_dir=f"{args.output_dir}/adapters")
        adapter_files = []
        for dup_file in duplicate_files:
            result = adapter_remover.remove_adapters(dup_file)
            if result:
                adapter_files.append(result)
        logger.info(f"Completed adapter removal. Generated {len(adapter_files)} files")
        
        # 3. Quality Filtering
        logger.info("\n=== Step 3: Quality Filtering ===")
        quality_filter = QualityFilter(output_dir=f"{args.output_dir}/quality")
        quality_files = []
        for adapter_file in adapter_files:
            result = quality_filter.filter_sequences(adapter_file)
            if result:
                quality_files.append(result)
        logger.info(f"Completed quality filtering. Generated {len(quality_files)} files")
        
        # 4. Pair Matching
        logger.info("\n=== Step 4: Pair Matching ===")
        pair_matcher = PairMatcher(output_dir=f"{args.output_dir}/paired", num_cores=num_cores)
        paired_files = pair_matcher.process_pairs(quality_files)
        logger.info(f"Completed pair matching. Generated {len(paired_files)} paired files")
        
        # 5. Secondary Duplicate Removal
        logger.info("\n=== Step 5: Secondary Duplicate Removal ===")
        secondary_duplicate_remover = DuplicateRemover(output_dir=f"{args.output_dir}/new_duplicate")
        final_files = []
        for paired_file in paired_files:
            result = secondary_duplicate_remover.remove_duplicates(paired_file)
            if result:
                final_files.append(result)
        logger.info(f"Completed secondary duplicate removal. Generated {len(final_files)} files")
        
        # 6. Alignment Analysis
        logger.info("\n=== Step 6: Alignment Analysis ===")
        alignment_analyzer = AlignmentAnalyzer(
            input_dir=f"{args.output_dir}/new_duplicate",
            output_dir=f"{args.output_dir}/alignment",
            num_cores=num_cores
        )
        alignment_files = alignment_analyzer.process_files()
        logger.info(f"Completed alignment analysis. Generated {len(alignment_files)} files")
        
        # 7. BLAST Search and Report Generation
        logger.info("\n=== Step 7: BLAST Search ===")
        blast_searcher = BlastSearcher(
            input_dir=f"{args.output_dir}/alignment",
            output_dir=f"{args.output_dir}/blast"
        )
        blast_files = blast_searcher.process_files()
        
        if blast_files:
            logger.info("Generating BLAST PDF report...")
            pdf_generator = BlastPDFGenerator(
                blast_dir=f"{args.output_dir}/blast",
                output_dir=f"{args.output_dir}/reports"
            )
            pdf_report = pdf_generator.generate_report()
            if pdf_report:
                logger.info(f"PDF report generated: {pdf_report}")
        
        # Print final summary
        logger.info("\n=== Pipeline Summary ===")
        logger.info(f"Input files: {len(input_files)}")
        logger.info(f"After initial duplicate removal: {len(duplicate_files)}")
        logger.info(f"After adapter removal: {len(adapter_files)}")
        logger.info(f"After quality filtering: {len(quality_files)}")
        logger.info(f"After pair matching: {len(paired_files)}")
        logger.info(f"After secondary duplicate removal: {len(final_files)}")
        logger.info(f"Alignment files generated: {len(alignment_files)}")
        logger.info(f"BLAST results generated: {len(blast_files) if blast_files else 0}")
        
        logger.info("\nPipeline completed successfully!")
        
    except Exception as e:
        logger.error(f"Pipeline error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main() 