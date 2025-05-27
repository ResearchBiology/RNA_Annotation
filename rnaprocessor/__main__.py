import argparse
import logging
import shutil
from pathlib import Path
from new_rnaprocessor.src.sequence_cleaner import SequenceCleaner
from new_rnaprocessor.src.adapter_remover import AdapterRemover
from new_rnaprocessor.src.quality_filter import QualityFilter
from new_rnaprocessor.src.pair_filter import PairFilter
from new_rnaprocessor.src.alignment_analyzer import AlignmentAnalyzer
from new_rnaprocessor.src.blast_searcher import BlastSearcher

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def setup_directories(base_output_dir):
    """Setup output directories and return paths"""
    base_dir = Path(base_output_dir)
    
    # Only create directories we want to keep outputs from
    final_dirs = {
        'cleaned': base_dir / 'final_cleaned',
        'aligned': base_dir / 'alignment_ready',
        'blast': base_dir / 'blast_results'
    }
    
    # Create final output directories
    for dir_path in final_dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)
    
    # Create temp directory for intermediate files
    temp_dir = base_dir / 'temp'
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    return final_dirs, temp_dir

def cleanup_temp_files(temp_dir):
    """Remove temporary directory and its contents"""
    try:
        shutil.rmtree(temp_dir)
        logger.info("Cleaned up temporary files")
    except Exception as e:
        logger.warning(f"Error cleaning up temporary files: {str(e)}")

def process_files(input_dir, output_dir, email=None):
    """Process all input files through the RNA processing pipeline"""
    try:
        # Setup directories
        final_dirs, temp_dir = setup_directories(output_dir)
        
        # Initialize processors
        cleaner = SequenceCleaner(output_dir=str(temp_dir / "cleaned"))
        adapter_remover = AdapterRemover(output_dir=str(temp_dir / "adapter_removed"))
        quality_filter = QualityFilter(output_dir=str(temp_dir / "quality_filtered"))
        pair_filter = PairFilter(output_dir=str(temp_dir / "paired"))
        final_cleaner = SequenceCleaner(output_dir=str(final_dirs['cleaned']))
        alignment_analyzer = AlignmentAnalyzer(output_dir=str(final_dirs['aligned']))
        
        if email:
            blast_searcher = BlastSearcher(
                output_dir=str(final_dirs['blast']),
                email=email
            )
        
        # Process each pair of input files
        input_path = Path(input_dir)
        forward_files = sorted(input_path.glob("*_1.fq.gz"))
        
        for forward_file in forward_files:
            base_name = forward_file.name[:-8]  # Remove _1.fq.gz
            reverse_file = input_path / f"{base_name}_2.fq.gz"
            
            if not reverse_file.exists():
                logger.warning(f"No matching reverse file for {forward_file.name}")
                continue
            
            logger.info(f"\nProcessing pair: {forward_file.name} and {reverse_file.name}")
            
            try:
                # Step 1: Initial cleaning
                logger.info("\n=== Step 1: Initial Sequence Cleaning ===")
                cleaned_forward = cleaner.process_file(str(forward_file))
                cleaned_reverse = cleaner.process_file(str(reverse_file))
                
                # Step 2: Adapter removal
                logger.info("\n=== Step 2: Adapter Removal ===")
                adapter_removed_forward = adapter_remover.remove_adapters(cleaned_forward)
                adapter_removed_reverse = adapter_remover.remove_adapters(cleaned_reverse)
                
                # Step 3: Quality filtering
                logger.info("\n=== Step 3: Quality Filtering ===")
                quality_filtered_forward = quality_filter.filter_sequences(adapter_removed_forward)
                quality_filtered_reverse = quality_filter.filter_sequences(adapter_removed_reverse)
                
                # Step 4: Pair filtering
                logger.info("\n=== Step 4: Pair Filtering ===")
                paired_file, _ = pair_filter.combine_paired_files(
                    quality_filtered_forward,
                    quality_filtered_reverse
                )
                
                # Step 5: Final cleaning (duplicate removal)
                logger.info("\n=== Step 5: Final Cleaning ===")
                final_cleaned = final_cleaner.process_file(paired_file)
                
                # Step 6: Alignment analysis
                logger.info("\n=== Step 6: Alignment Analysis ===")
                fasta_file, _ = alignment_analyzer.analyze_file(final_cleaned)
                
                # Step 7: BLAST search (if email provided)
                if email and fasta_file:
                    logger.info("\n=== Step 7: BLAST Search ===")
                    blast_searcher.search(fasta_file)
                
            except Exception as e:
                logger.error(f"Error processing {base_name}: {str(e)}")
                continue
        
        # Cleanup temporary files
        cleanup_temp_files(temp_dir)
        
        logger.info("\nPipeline completed successfully!")
        logger.info(f"Final outputs can be found in:")
        logger.info(f"- Cleaned sequences: {final_dirs['cleaned']}")
        logger.info(f"- Alignment results: {final_dirs['aligned']}")
        if email:
            logger.info(f"- BLAST results: {final_dirs['blast']}")
        
    except Exception as e:
        logger.error(f"Pipeline error: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description="RNA Processing Pipeline")
    
    parser.add_argument(
        '--input-dir',
        type=str,
        required=True,
        help='Directory containing input .fq.gz files'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Directory for output files'
    )
    
    parser.add_argument(
        '--email',
        type=str,
        help='Email for NCBI BLAST (optional)'
    )
    
    args = parser.parse_args()
    
    # Run the pipeline
    process_files(args.input_dir, args.output_dir, args.email)

if __name__ == "__main__":
    main()
