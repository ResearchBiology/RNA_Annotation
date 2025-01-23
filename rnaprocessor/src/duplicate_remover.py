import logging
from pathlib import Path

logger = logging.getLogger(__name__)

class DuplicateRemover:
    """Class for removing duplicate sequences from paired FASTQ files"""
    
    def __init__(self, output_dir="results/new_duplicate"):
        """
        Initialize DuplicateRemover
        
        Args:
            output_dir (str): Directory to store unique sequences
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def remove_duplicates(self, input_file):
        """
        Remove duplicate sequence pairs from paired .fq file
        
        Args:
            input_file (str): Path to input paired .fq file
            
        Returns:
            str: Path to processed file
        """
        try:
            # Setup paths
            input_path = Path(input_file)
            output_path = self.output_dir / f"unique_{input_path.name}"
            stats_path = self.output_dir / f"{input_path.stem}_stats.txt"
            
            # Track unique sequence pairs
            seen_pairs = set()
            total_pairs = 0
            unique_pairs = 0
            
            logger.info(f"Processing file: {input_path.name}")
            
            # Process file
            with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
                while True:
                    # Read four lines (two sequence pairs)
                    info1 = infile.readline().strip()
                    seq1 = infile.readline().strip()
                    info2 = infile.readline().strip()
                    seq2 = infile.readline().strip()
                    
                    # Check for end of file
                    if not info1 or not seq1 or not info2 or not seq2:
                        break
                        
                    total_pairs += 1
                    
                    # Create a tuple of both sequences for uniqueness check
                    pair_key = (seq1, seq2)
                    
                    # If pair is unique, write both sequences
                    if pair_key not in seen_pairs:
                        seen_pairs.add(pair_key)
                        unique_pairs += 1
                        outfile.write(f"{info1}\n{seq1}\n{info2}\n{seq2}\n")
                    
                    # Log progress
                    if total_pairs % 10000 == 0:
                        logger.info(f"Processed {total_pairs} sequence pairs")
            
            # Calculate statistics
            duplicates = total_pairs - unique_pairs
            
            # Save statistics
            with open(stats_path, 'w') as f:
                f.write("=== Duplicate Removal Statistics ===\n\n")
                f.write(f"Input file: {input_path.name}\n")
                f.write(f"Total sequence pairs: {total_pairs}\n")
                f.write(f"Unique sequence pairs: {unique_pairs}\n")
                f.write(f"Duplicate pairs removed: {duplicates}\n")
                f.write(f"Duplication rate: {(duplicates/total_pairs)*100:.2f}%\n")
            
            logger.info(f"\nProcessing completed:")
            logger.info(f"- Total pairs: {total_pairs}")
            logger.info(f"- Unique pairs: {unique_pairs}")
            logger.info(f"- Duplicates removed: {duplicates}")
            logger.info(f"- Results saved to: {output_path}")
            logger.info(f"- Statistics saved to: {stats_path}")
            
            return str(output_path)
            
        except Exception as e:
            logger.error(f"Error processing file {input_file}: {str(e)}")
            raise 