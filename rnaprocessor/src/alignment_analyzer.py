import logging
from pathlib import Path
from Bio import Align
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from collections import Counter

logger = logging.getLogger(__name__)

class AlignmentAnalyzer:
    """Class for analyzing sequence alignments"""
    
    def __init__(self, input_dir="results/new_duplicate", output_dir="results/alignment", num_cores=22):
        """Initialize AlignmentAnalyzer"""
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.num_cores = min(num_cores, multiprocessing.cpu_count())
    
    def _read_sequences(self, file_path):
        """Read sequences from file, returning list of sequences"""
        sequences = []
        with open(file_path, 'r') as f:
            while True:
                info_line = f.readline().strip()
                seq_line = f.readline().strip()
                
                if not info_line or not seq_line:
                    break
                    
                sequences.append((info_line, seq_line))
        return sequences
    
    def _process_file(self, file_path):
        """Process sequences in a single file and save to txt"""
        try:
            logger.info(f"Processing {file_path.name}")
            
            # Get sample name
            sample_name = file_path.stem.replace("unique_paired_quality_filtered_processed_unique_cleaned_", "")
            sample_name = sample_name.replace("_L2", "")
            
            # Read sequences
            sequences = self._read_sequences(file_path)
            
            # Find most common sequence
            if sequences:
                # Count sequence occurrences
                seq_counter = Counter(seq for _, seq in sequences)
                most_common_seq = seq_counter.most_common(1)[0]  # Returns (sequence, count)
                
                # Find the first info line for this sequence
                info_line = next(info for info, seq in sequences if seq == most_common_seq[0])
                
                output_file = self.output_dir / f"{sample_name}_sequence.txt"
                
                # Write to file
                with open(output_file, 'w') as f:
                    f.write(f"Sample: {sample_name}\n")
                    f.write(f"Sequence ID: {info_line}\n")
                    f.write(f"Length: {len(most_common_seq[0])}\n")
                    f.write(f"Frequency: {most_common_seq[1]} times\n")
                    f.write(f"Sequence:\n{most_common_seq[0]}\n")
                
                logger.info(f"Saved most common sequence (length {len(most_common_seq[0])}, "
                          f"frequency {most_common_seq[1]}) for {sample_name}")
                return output_file
            else:
                logger.warning(f"No sequences found in {file_path.name}")
                return None
                
        except Exception as e:
            logger.error(f"Error processing {file_path.name}: {str(e)}")
            return None
    
    def process_files(self):
        """Process all files and save individual txt files"""
        # Get list of files
        input_files = list(self.input_dir.glob("unique_paired_*.fq"))
        processed_files = []
        
        # Process files using ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=self.num_cores) as executor:
            futures = [executor.submit(self._process_file, file_path) for file_path in input_files]
            for future in futures:
                result = future.result()
                if result:
                    processed_files.append(result)
        
        logger.info(f"\nProcessing completed:")
        logger.info(f"- Total files processed: {len(processed_files)}")
        logger.info(f"- Results saved in: {self.output_dir}")
        
        return processed_files 