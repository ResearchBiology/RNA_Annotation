import logging
import gzip
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class SequenceCleaner:
    """Class for cleaning FASTQ.GZ files and removing duplicates while maintaining FASTQ format"""
    
    def __init__(self, output_dir="cleaned_sequences"):
        """
        Initialize SequenceCleaner
        
        Args:
            output_dir (str): Directory to store cleaned sequences
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def _calculate_sequence_hash(self, sequence, quality):
        """Calculate a hash that considers both sequence and quality"""
        return hash(sequence + '|' + quality)

    def process_file(self, input_file):
        """
        Process a FASTQ file and remove duplicates
        
        Args:
            input_file (str): Path to input .fq file or .fq.gz file
            
        Returns:
            str: Path to cleaned FASTQ file
        """
        input_path = Path(input_file)
        output_path = self.output_dir / f"{input_path.stem}_unique.fq"
        stats_path = self.output_dir / f"{input_path.stem}_unique_stats.txt"
        
        try:
            logger.info(f"Processing file: {input_path.name}")
            
            # Track unique sequences
            seen_hashes = {}
            total_sequences = 0
            duplicate_sequences = 0
            
            # Process sequences
            logger.info("Reading sequences and removing duplicates...")
            
            # Determine if the file is gzipped
            is_gzipped = input_path.suffix.lower() == '.gz'
            
            # Open file with appropriate method
            opener = gzip.open if is_gzipped else open
            
            with opener(input_path, 'rt') as f:  # 'rt' mode for text reading
                # Process 4 lines at a time (FASTQ format)
                while True:
                    # Read a complete FASTQ record
                    header = f.readline().strip()
                    if not header:  # End of file
                        break
                        
                    sequence = f.readline().strip()
                    plus_line = f.readline().strip()
                    quality = f.readline().strip()
                    
                    total_sequences += 1
                    
                    # Calculate hash considering both sequence and quality
                    seq_hash = self._calculate_sequence_hash(sequence, quality)
                    
                    # Only keep first occurrence of each sequence
                    if seq_hash not in seen_hashes:
                        seen_hashes[seq_hash] = (header, sequence, plus_line, quality)
                    else:
                        duplicate_sequences += 1
                    
            
            # Write unique sequences to FASTQ
            logger.info("Writing unique sequences...")
            with open(output_path, 'w') as f:
                for header, sequence, plus_line, quality in seen_hashes.values():
                    f.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")
            
            # Calculate statistics
            unique_sequences = len(seen_hashes)
            duplicate_rate = (duplicate_sequences) / total_sequences * 100
            
            # Save statistics
            logger.info("Saving statistics...")
            with open(stats_path, 'w') as f:
                f.write("=== Final Sequence Cleaning Statistics ===\n\n")
                f.write(f"Input file: {input_path.name}\n\n")
                f.write("Summary:\n")
                f.write(f"Total sequences processed: {total_sequences}\n")
                f.write(f"Unique sequences kept: {unique_sequences}\n")
                f.write(f"Duplicate sequences removed: {duplicate_sequences}\n")
                f.write(f"Duplication rate: {duplicate_rate:.2f}%\n")
            
            logger.info(f"\nCleaning completed:")
            logger.info(f"- Total sequences: {total_sequences}")
            logger.info(f"- Unique sequences: {unique_sequences}")
            logger.info(f"- Duplicates removed: {duplicate_sequences}")
            logger.info(f"- Results saved to: {output_path}")
            logger.info(f"- Statistics saved to: {stats_path}")
            
            return str(output_path)
            
        except Exception as e:
            logger.error(f"Error processing file {input_file}: {str(e)}")
            raise
            
    def process_directory(self, input_dir):
        """
        Process all FASTQ files in a directory (.fq, .fastq, .fq.gz, .fastq.gz)
        
        Args:
            input_dir (str): Directory containing FASTQ files
            
        Returns:
            list: Paths to processed files
        """
        input_path = Path(input_dir)
        processed_files = []
        
        try:
            # Get all FASTQ files (both gzipped and non-gzipped)
            fastq_files = []
            fastq_files.extend(list(input_path.glob("*.fq.gz")))
            fastq_files.extend(list(input_path.glob("*.fastq.gz")))
            fastq_files.extend(list(input_path.glob("*.fq")))
            fastq_files.extend(list(input_path.glob("*.fastq")))
            
            logger.info(f"Found {len(fastq_files)} FASTQ files to process")
            
            for fastq_file in fastq_files:
                try:
                    processed_file = self.process_file(fastq_file)
                    processed_files.append(processed_file)
                except Exception as e:
                    logger.error(f"Error processing {fastq_file.name}: {str(e)}")
                    continue
                    
            return processed_files
            
        except Exception as e:
            logger.error(f"Error processing directory {input_dir}: {str(e)}")
            raise
