import logging
from pathlib import Path
from Bio.Seq import Seq
from collections import defaultdict
from tqdm import tqdm
import concurrent.futures

logger = logging.getLogger(__name__)

class PairFilter:
    """Class for pairing and combining forward and reverse sequences"""
    
    def __init__(self, output_dir="paired_sequences", num_cores=6):
        """
        Initialize PairFilter
        
        Args:
            output_dir (str): Directory to store paired sequences
            num_cores (int): Number of CPU cores to use
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.num_cores = num_cores
    
    def _read_sequences(self, file_path):
        """
        Read all sequences from a file
        
        Args:
            file_path (str): Path to the sequence file
            
        Returns:
            dict: Dictionary of sequences {id: (header, sequence, plus_line, quality)}
        """
        sequences = {}
        with open(file_path) as f:
            while True:
                header = f.readline().rstrip('\n')
                if not header:
                    break
                sequence = f.readline().rstrip('\n')
                plus_line = f.readline().rstrip('\n')
                quality = f.readline().rstrip('\n')
                
                if header and sequence:
                    # Extract sequence ID from Illumina format
                    # Format: @INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y READ:FILTERED:FLAG:INDEX
                    seq_id = header[1:].split()[0]  # Get everything before the space
                    # Remove the read number (1 or 2) from comparison
                    seq_id = seq_id.rsplit(' ', 1)[0] if ' ' in seq_id else seq_id
                    sequences[seq_id] = (header, sequence, plus_line, quality)
        
        return sequences
    
    def _calculate_similarity(self, seq1, seq2):
        """Calculate similarity between two sequences"""
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / max(len(seq1), len(seq2))

    def _is_complementary(self, seq1, seq2):
        """Check if sequences are complementary"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        rev_comp = ''.join(complement.get(base, 'N') for base in reversed(seq2))
        similarity = self._calculate_similarity(seq1, rev_comp)
        return similarity >= 0.8

    def _find_best_match(self, forward_seq, reverse_seqs):
        """Find best matching sequence based on complementarity or similarity"""
        best_match = None
        best_score = 0

        for r_id, (r_header, r_seq, r_plus, r_qual) in reverse_seqs.items():
            # Try complementary matching first
            if self._is_complementary(forward_seq, r_seq):
                return r_id, (r_header, r_seq, r_plus, r_qual)
            
            # If not complementary, check similarity
            similarity = self._calculate_similarity(forward_seq, r_seq)
            if similarity >= 0.8 and similarity > best_score:
                best_score = similarity
                best_match = (r_id, (r_header, r_seq, r_plus, r_qual))

        return best_match

    def _process_chunk(self, chunk_data):
        """Process a chunk of sequences using sequence content matching"""
        forward_chunk, reverse_seqs = chunk_data
        pairs = []
        unpaired = set()
        used_reverse = set()

        for f_id, (f_header, f_seq, f_plus, f_qual) in forward_chunk.items():
            # Find best matching reverse sequence
            match = self._find_best_match(f_seq, {k: v for k, v in reverse_seqs.items() if k not in used_reverse})
            
            if match:
                r_id, (r_header, r_seq, r_plus, r_qual) = match
                pairs.append((
                    f_header, f_seq, f_plus, f_qual,
                    r_header, r_seq, r_plus, r_qual
                ))
                used_reverse.add(r_id)
            else:
                unpaired.add(f_id)

        return pairs, unpaired
    
    def combine_paired_files(self, forward_file, reverse_file):
        """
        Combine paired sequences from forward and reverse files
        
        Args:
            forward_file (str): Path to forward reads file (_1.fq)
            reverse_file (str): Path to reverse reads file (_2.fq)
            
        Returns:
            tuple: (path to paired sequences file, path to stats file)
        """
        forward_path = Path(forward_file)
        reverse_path = Path(reverse_file)
        
        # Create output files
        base_name = forward_path.stem.replace('_1', '')
        paired_path = self.output_dir / f"{base_name}.fq"
        stats_path = self.output_dir / f"{base_name}_pairing_stats.txt"
        
        try:
            logger.info(f"Processing files: {forward_path.name} and {reverse_path.name}")
            
            # Read all sequences
            logger.info("Reading sequences...")
            forward_seqs = self._read_sequences(forward_file)
            reverse_seqs = self._read_sequences(reverse_file)
            
            # Process in chunks for memory efficiency
            chunk_size = max(1000, len(forward_seqs) // (self.num_cores * 10))
            forward_chunks = [
                dict(list(forward_seqs.items())[i:i + chunk_size])
                for i in range(0, len(forward_seqs), chunk_size)
            ]
            
            # Process chunks in parallel
            pairs = []
            unpaired_forward = set()
            
            with concurrent.futures.ProcessPoolExecutor(max_workers=self.num_cores) as executor:
                futures = []
                for chunk in forward_chunks:
                    futures.append(
                        executor.submit(self._process_chunk, (chunk, reverse_seqs))
                    )
                
                for future in tqdm(concurrent.futures.as_completed(futures),
                                 total=len(futures),
                                 desc="Processing chunks"):
                    chunk_pairs, chunk_unpaired = future.result()
                    pairs.extend(chunk_pairs)
                    unpaired_forward.update(chunk_unpaired)
            
            # Write paired sequences
            with open(paired_path, 'w') as f:
                for (f_header, f_seq, f_plus, f_qual,
                     r_header, r_seq, r_plus, r_qual) in pairs:
                    # Write forward read
                    f.write(f"{f_header}\n{f_seq}\n{f_plus}\n{f_qual}\n")
                    # Write reverse read
                    f.write(f"{r_header}\n{r_seq}\n{r_plus}\n{r_qual}\n")
            
            # Calculate statistics
            total_forward = len(forward_seqs)
            total_reverse = len(reverse_seqs)
            total_pairs = len(pairs)
            unpaired_reverse = set(reverse_seqs.keys()) - set(forward_seqs.keys())
            
            # Save statistics
            with open(stats_path, 'w') as f:
                f.write("=== Sequence Pairing Statistics ===\n\n")
                f.write(f"Forward file: {forward_path.name}\n")
                f.write(f"Reverse file: {reverse_path.name}\n\n")
                f.write("Summary:\n")
                f.write(f"Total forward sequences: {total_forward}\n")
                f.write(f"Total reverse sequences: {total_reverse}\n")
                f.write(f"Paired sequences found: {total_pairs}\n")
                f.write(f"Pairing rate: {total_pairs/max(total_forward, total_reverse)*100:.2f}%\n")
                
                if unpaired_forward:
                    f.write("\nFirst 10 sequences found only in forward file:\n")
                    for seq_id in sorted(list(unpaired_forward))[:10]:
                        f.write(f"- {seq_id}\n")
                    if len(unpaired_forward) > 10:
                        f.write(f"... and {len(unpaired_forward) - 10} more\n")
                
                if unpaired_reverse:
                    f.write("\nFirst 10 sequences found only in reverse file:\n")
                    for seq_id in sorted(list(unpaired_reverse))[:10]:
                        f.write(f"- {seq_id}\n")
                    if len(unpaired_reverse) > 10:
                        f.write(f"... and {len(unpaired_reverse) - 10} more\n")
            
            logger.info(f"\nPairing completed:")
            logger.info(f"- Total forward sequences: {total_forward}")
            logger.info(f"- Total reverse sequences: {total_reverse}")
            logger.info(f"- Paired sequences: {total_pairs}")
            logger.info(f"- Results saved to: {paired_path}")
            logger.info(f"- Statistics saved to: {stats_path}")
            
            return str(paired_path), str(stats_path)
            
        except Exception as e:
            logger.error(f"Error during sequence pairing: {str(e)}")
            raise
