import logging
from pathlib import Path
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

logger = logging.getLogger(__name__)

class PairMatcher:
    """Match and process forward and reverse read pairs"""
    
    def __init__(self, output_dir="results/paired", num_cores=22):
        """Initialize PairMatcher"""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.reversed_dir = self.output_dir / "reversed"
        self.reversed_dir.mkdir(exist_ok=True)
        self.min_overlap = 30  # Minimum overlap length 1/5
        self.num_cores = min(num_cores, multiprocessing.cpu_count())
    
    def _find_sequence_overlap(self, seq1, seq2):
        """
        Check if seq2 contains part of seq1 or vice versa.
        Looking for significant overlapping regions.
        """
        # Try to find seq1 within seq2
        for i in range(len(seq1) - self.min_overlap + 1):
            for j in range(i + self.min_overlap, len(seq1) + 1):
                substring = seq1[i:j]
                if substring in seq2:
                    return True
        return False
    
    def _process_chunk(self, chunk_data):
        """Process a chunk of sequences"""
        forward_sequences, reversed_data = chunk_data
        matched_pairs = []
        
        for f_info, f_seq in forward_sequences:
            # Look for overlapping sequences
            for rev_info, rev_seq in reversed_data:
                if self._find_sequence_overlap(f_seq, rev_seq):
                    matched_pairs.append((f_info, f_seq, rev_info, rev_seq))
                    break
        
        return matched_pairs
    
    def _identify_pairs(self, file_list):
        """Identify matching forward and reverse read files"""
        pairs = []
        forward_files = {}
        
        for file_path in file_list:
            path = Path(file_path)
            if '_1.fq' in str(path):
                base_name = str(path).replace('_1.fq', '')
                forward_files[base_name] = path
        
        for file_path in file_list:
            path = Path(file_path)
            if '_2.fq' in str(path):
                base_name = str(path).replace('_2.fq', '')
                if base_name in forward_files:
                    pairs.append((forward_files[base_name], path))
        
        return pairs
    
    def _reverse_sequence(self, sequence):
        """Simply reverse the sequence"""
        return sequence[::-1]
    
    def process_pairs(self, input_files):
        """Process and pair forward and reverse reads"""
        try:
            pairs = self._identify_pairs(input_files)
            if not pairs:
                raise ValueError("No valid pairs found in input files")
            
            processed_files = []
            
            for forward_file, reverse_file in pairs:
                logger.info(f"\nProcessing pair:")
                logger.info(f"Forward: {forward_file.name}")
                logger.info(f"Reverse: {reverse_file.name}")
                
                # Setup output paths
                base_name = forward_file.stem.replace('_1', '')
                output_path = self.output_dir / f"paired_{base_name}.fq"
                reversed_path = self.reversed_dir / f"reversed_{reverse_file.name}"
                stats_path = self.output_dir / f"{base_name}_pairing_stats.txt"
                
                # Read sequences from both files
                forward_sequences = []
                with open(forward_file, 'r') as f:
                    while True:
                        info_line = f.readline().strip()
                        seq_line = f.readline().strip()
                        if not info_line or not seq_line:
                            break
                        forward_sequences.append((info_line, seq_line))
                
                reversed_data = []
                with open(reverse_file, 'r') as f:
                    while True:
                        info_line = f.readline().strip()
                        seq_line = f.readline().strip()
                        if not info_line or not seq_line:
                            break
                        reversed_data.append((info_line, seq_line))
                
                # Split into chunks for parallel processing
                chunk_size = max(1, len(forward_sequences) // self.num_cores)
                chunks = [
                    (forward_sequences[i:i + chunk_size], reversed_data)
                    for i in range(0, len(forward_sequences), chunk_size)
                ]
                
                # Process chunks in parallel
                matched_pairs = []
                with ProcessPoolExecutor(max_workers=self.num_cores) as executor:
                    futures = [executor.submit(self._process_chunk, chunk) for chunk in chunks]
                    for future in futures:
                        matched_pairs.extend(future.result())
                
                # Write results
                with open(output_path, 'w') as out_file:
                    for f_info, f_seq, r_info, r_seq in matched_pairs:
                        out_file.write(f"{f_info}\n{f_seq}\n")
                        out_file.write(f"{r_info}\n{r_seq}\n")
                
                # Calculate statistics
                total_sequences = len(forward_sequences)
                paired_sequences = len(matched_pairs)
                
                # Save statistics
                with open(stats_path, 'w') as f:
                    f.write("=== Sequence Pairing Statistics ===\n\n")
                    f.write(f"Forward file: {forward_file.name}\n")
                    f.write(f"Reverse file: {reverse_file.name}\n\n")
                    f.write("Summary:\n")
                    f.write(f"Total sequences processed: {total_sequences}\n")
                    f.write(f"Sequences paired: {paired_sequences}\n")
                    f.write(f"Pairing rate: {paired_sequences/total_sequences*100:.2f}%\n")
                    f.write("\nPairing criteria:\n")
                    f.write(f"- Minimum overlap length: {self.min_overlap} bases\n")
                    f.write(f"\nProcessing details:\n")
                    f.write(f"- Cores used: {self.num_cores}\n")
                    f.write(f"- Chunk size: {chunk_size} sequences\n")
                
                logger.info(f"\nPairing completed:")
                logger.info(f"- Total sequences: {total_sequences}")
                logger.info(f"- Paired sequences: {paired_sequences}")
                logger.info(f"- Using {self.num_cores} cores")
                logger.info(f"- Results saved to: {output_path}")
                logger.info(f"- Statistics saved to: {stats_path}")
                
                processed_files.append(str(output_path))
            
            return processed_files
            
        except Exception as e:
            logger.error(f"Error during pair processing: {str(e)}")
            raise 