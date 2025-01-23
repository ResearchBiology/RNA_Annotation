import logging
from pathlib import Path

logger = logging.getLogger(__name__)

class AdapterRemover:
    """Class for removing adapter sequences based on identifier tag"""
    
    def __init__(self, output_dir="results/adapters"):
        """
        Initialize AdapterRemover
        
        Args:
            output_dir (str): Directory to store processed files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.identifier = "ATAGCGG"
        self.check_region = 30  # Last 30bp to check for identifier
    
    def remove_adapters(self, input_file):
        """
        Remove sequences based on identifier presence in last 30bp.
        If identifier found, trim sequence up to identifier.
        If not found in last 30bp, remove entire sequence.
        
        Args:
            input_file (str): Path to input file
            
        Returns:
            str: Path to processed file
        """
        try:
            # Setup paths
            input_path = Path(input_file)
            output_path = self.output_dir / f"processed_{input_path.name}"
            stats_path = self.output_dir / f"{input_path.stem}_adapter_stats.txt"
            
            # Track statistics
            total_sequences = 0
            trimmed_sequences = 0
            removed_sequences = 0
            
            logger.info(f"Processing file: {input_path.name}")
            
            # Process file
            with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
                while True:
                    # Read pair of lines
                    info_line = infile.readline().strip()
                    seq_line = infile.readline().strip()
                    
                    # Check for end of file
                    if not info_line or not seq_line:
                        break
                        
                    total_sequences += 1
                    
                    # Check last 30bp for identifier
                    check_start = max(0, len(seq_line) - self.check_region)
                    last_region = seq_line[check_start:]
                    
                    # Find identifier position in last region
                    identifier_pos = last_region.find(self.identifier)
                    
                    if identifier_pos != -1:
                        # Identifier found - trim sequence
                        actual_pos = check_start + identifier_pos
                        trimmed_seq = seq_line[:actual_pos]
                        outfile.write(f"{info_line}\n{trimmed_seq}\n")
                        trimmed_sequences += 1
                    else:
                        # Identifier not found in last 30bp - remove sequence
                        removed_sequences += 1
                    
                    # Log progress
                    if total_sequences % 100000 == 0:
                        logger.info(f"Processed {total_sequences} sequences")
            
            # Save statistics
            with open(stats_path, 'w') as f:
                f.write("=== Adapter Removal Statistics ===\n\n")
                f.write(f"Input file: {input_path.name}\n")
                f.write(f"Identifier tag: {self.identifier}\n")
                f.write(f"Check region: last {self.check_region} bp\n\n")
                f.write("Summary:\n")
                f.write(f"Total sequences: {total_sequences}\n")
                f.write(f"Sequences trimmed: {trimmed_sequences} ({trimmed_sequences/total_sequences*100:.2f}%)\n")
                f.write(f"Sequences removed: {removed_sequences} ({removed_sequences/total_sequences*100:.2f}%)\n")
                f.write(f"Sequences kept: {trimmed_sequences} ({trimmed_sequences/total_sequences*100:.2f}%)\n")
            
            logger.info(f"\nProcessing completed:")
            logger.info(f"- Total sequences: {total_sequences}")
            logger.info(f"- Sequences trimmed: {trimmed_sequences}")
            logger.info(f"- Sequences removed: {removed_sequences}")
            logger.info(f"- Results saved to: {output_path}")
            logger.info(f"- Statistics saved to: {stats_path}")
            
            return str(output_path)
            
        except Exception as e:
            logger.error(f"Error processing file {input_file}: {str(e)}")
            raise 