import logging
from pathlib import Path
import statistics
import concurrent.futures
from tqdm import tqdm

logger = logging.getLogger(__name__)

class QualityFilter:
    """Class for filtering sequences based on quality scores"""
    
    def __init__(self, output_dir="quality_filtered", min_quality=15, window_size=30, num_cores=6):
        """
        Initialize QualityFilter
        
        Args:
            output_dir (str): Directory to store filtered sequences
            min_quality (int): Minimum average quality score
            window_size (int): Window size for quality calculation
            num_cores (int): Number of CPU cores to use
        """
        self.output_dir = Path(output_dir)
        self.min_quality = min_quality
        self.window_size = window_size
        self.num_cores = num_cores
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def _check_quality(self, quality_str):
        """
        Check if sequence passes quality threshold using sliding window
        
        Args:
            quality_str (str): Quality string in ASCII format
            
        Returns:
            bool: True if sequence passes quality check, False otherwise
        """
        # Convert ASCII to quality scores (Phred+33)
        quality_scores = [ord(c) - 33 for c in quality_str]
        
        # Calculate average quality in sliding windows
        for i in range(0, len(quality_scores) - self.window_size + 1):
            window = quality_scores[i:i + self.window_size]
            if statistics.mean(window) < self.min_quality:
                return False
        return True
    
    def _process_chunk(self, chunk_data):
        """Process a chunk of sequences"""
        filtered_sequences = []
        stats = {"total": 0, "passed": 0}
        
        for header, sequence, plus_line, quality in chunk_data:
            stats["total"] += 1
            if self._check_quality(quality):
                filtered_sequences.append((header, sequence, plus_line, quality))
                stats["passed"] += 1
        
        return filtered_sequences, stats
    
    def filter_sequences(self, input_file):
        """
        Filter sequences based on quality scores using parallel processing
        
        Args:
            input_file (str): Path to input file
            
        Returns:
            str: Path to filtered file
        """
        input_path = Path(input_file)
        output_path = self.output_dir / f"{input_path.name}"
        stats_path = self.output_dir / f"{input_path.stem}_quality_stats.txt"
        
        try:
            logger.info(f"Processing file: {input_path.name}")
            
            # Read sequences in chunks
            sequences = []
            total_sequences = 0
            
            with open(input_file) as f:
                while True:
                    header = f.readline().rstrip('\n')
                    if not header:
                        break
                    sequence = f.readline().rstrip('\n')
                    plus_line = f.readline().rstrip('\n')
                    quality = f.readline().rstrip('\n')
                    
                    if header and sequence:
                        sequences.append((header, sequence, plus_line, quality))
                        total_sequences += 1
            
            # Process in parallel
            chunk_size = max(1000, total_sequences // (self.num_cores * 10))
            chunks = [sequences[i:i + chunk_size] for i in range(0, len(sequences), chunk_size)]
            
            filtered_sequences = []
            total_stats = {"total": 0, "passed": 0}
            
            with concurrent.futures.ProcessPoolExecutor(max_workers=self.num_cores) as executor:
                futures = []
                for chunk in chunks:
                    futures.append(executor.submit(self._process_chunk, chunk))
                
                for future in tqdm(concurrent.futures.as_completed(futures), 
                                 total=len(futures),
                                 desc="Processing chunks"):
                    chunk_sequences, chunk_stats = future.result()
                    filtered_sequences.extend(chunk_sequences)
                    # Update total stats
                    for key in total_stats:
                        total_stats[key] += chunk_stats[key]
            
            # Write filtered sequences
            with open(output_path, 'w') as f:
                for header, sequence, plus_line, quality in filtered_sequences:
                    f.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")
            
            # Save statistics
            with open(stats_path, 'w') as f:
                f.write("=== Quality Filter Statistics ===\n\n")
                f.write(f"Input file: {input_path.name}\n")
                f.write(f"Quality parameters:\n")
                f.write(f"- Minimum quality score: {self.min_quality}\n")
                f.write(f"- Window size: {self.window_size}\n\n")
                f.write("Summary:\n")
                f.write(f"Total sequences processed: {total_stats['total']}\n")
                f.write(f"Sequences passed: {total_stats['passed']} ({total_stats['passed']/total_stats['total']*100:.2f}%)\n")
                f.write(f"Sequences filtered: {total_stats['total'] - total_stats['passed']} ({(total_stats['total'] - total_stats['passed'])/total_stats['total']*100:.2f}%)\n")
            
            logger.info(f"\nQuality filtering completed:")
            logger.info(f"- Total sequences: {total_stats['total']}")
            logger.info(f"- Sequences passed: {total_stats['passed']}")
            logger.info(f"- Results saved to: {output_path}")
            logger.info(f"- Statistics saved to: {stats_path}")
            
            return str(output_path)
            
        except Exception as e:
            logger.error(f"Error processing file {input_file}: {str(e)}")
            raise
