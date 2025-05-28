import logging
from pathlib import Path
import concurrent.futures
from tqdm import tqdm

logger = logging.getLogger(__name__)

class AdapterRemover:
    """Class for removing adapter sequences based on identifier tag"""
    
    def __init__(self, output_dir="adapter_removed", num_cores=6):
        """
        Initialize AdapterRemover
        
        Args:
            output_dir (str): Directory to store processed files
            num_cores (int): Number of CPU cores to use
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.num_cores = num_cores
        self.identifier = "AGATCGG"
        self.check_region = 30  # Last 30bp to check for identifier
    
    def _process_sequence(self, header, sequence, plus_line, quality):
        """
        Process a single sequence
        
        Returns:
            tuple: (header, sequence, plus_line, quality) or None if sequence should be removed
        """
        seq_length = len(sequence)
        
        # Check last 30bp region
        check_start = max(0, seq_length - self.check_region)
        last_region = sequence[check_start:]
        
        # Find identifier in last region
        identifier_pos = last_region.find(self.identifier)
        if identifier_pos != -1:
            # Identifier found in last 30bp
            # Calculate actual position in full sequence
            actual_pos = check_start + identifier_pos
            # Keep sequence up to identifier
            return (header, sequence[:actual_pos], plus_line, quality[:actual_pos])
        else:
            # Identifier not found in last 30bp - remove sequence
            return None
    
    def _process_chunk(self, chunk_data):
        """Process a chunk of sequences"""
        processed_sequences = []
        stats = {"total": 0, "trimmed": 0, "removed": 0}
        
        for header, sequence, plus_line, quality in chunk_data:
            stats["total"] += 1
            result = self._process_sequence(header, sequence, plus_line, quality)
            
            if result:
                header, processed_seq, plus_line, processed_qual = result
                if len(processed_seq) < len(sequence):
                    stats["trimmed"] += 1
                processed_sequences.append((header, processed_seq, plus_line, processed_qual))
            else:
                stats["removed"] += 1
        
        return processed_sequences, stats
    
    def remove_adapters(self, input_file):
        """
        Remove adapter sequences from a file using parallel processing
        
        Args:
            input_file (str): Path to input file
            
        Returns:
            str: Path to processed file
        """
        input_path = Path(input_file)
        output_path = self.output_dir / f"{input_path.name}"
        stats_path = self.output_dir / f"{input_path.stem}_adapter_stats.txt"
        
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
            
            processed_sequences = []
            total_stats = {"total": 0, "trimmed": 0, "removed": 0}
            
            with concurrent.futures.ProcessPoolExecutor(max_workers=self.num_cores) as executor:
                futures = []
                for chunk in chunks:
                    futures.append(executor.submit(self._process_chunk, chunk))
                
                for future in tqdm(concurrent.futures.as_completed(futures), 
                                 total=len(futures),
                                 desc="Processing chunks"):
                    chunk_sequences, chunk_stats = future.result()
                    processed_sequences.extend(chunk_sequences)
                    # Update total stats
                    for key in total_stats:
                        total_stats[key] += chunk_stats[key]
            
            # Write processed sequences
            with open(output_path, 'w') as f:
                for header, sequence, plus_line, quality in processed_sequences:
                    f.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")
            
            # Save statistics
            with open(stats_path, 'w') as f:
                f.write("=== Adapter Removal Statistics ===\n\n")
                f.write(f"Input file: {input_path.name}\n")
                f.write(f"Identifier tag: {self.identifier}\n")
                f.write(f"Check region: last {self.check_region} bp\n\n")
                f.write("Summary:\n")
                f.write(f"Total sequences processed: {total_stats['total']}\n")
                f.write(f"Sequences trimmed: {total_stats['trimmed']} ({total_stats['trimmed']/total_stats['total']*100:.2f}%)\n")
                f.write(f"Sequences removed: {total_stats['removed']} ({total_stats['removed']/total_stats['total']*100:.2f}%)\n")
                f.write(f"Sequences kept: {len(processed_sequences)} ({len(processed_sequences)/total_stats['total']*100:.2f}%)\n")
            
            logger.info(f"\nAdapter removal completed:")
            logger.info(f"- Total sequences: {total_stats['total']}")
            logger.info(f"- Sequences trimmed: {total_stats['trimmed']}")
            logger.info(f"- Sequences removed: {total_stats['removed']}")
            logger.info(f"- Sequences kept: {len(processed_sequences)}")
            logger.info(f"- Results saved to: {output_path}")
            logger.info(f"- Statistics saved to: {stats_path}")
            
            return str(output_path)
            
        except Exception as e:
            logger.error(f"Error processing file {input_file}: {str(e)}")
            raise
