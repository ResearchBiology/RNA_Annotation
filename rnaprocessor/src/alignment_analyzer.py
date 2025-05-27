import logging
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
from tqdm import tqdm

logger = logging.getLogger(__name__)

class AlignmentAnalyzer:
    """Class for analyzing sequence alignments and preparing for BLAST"""
    
    def __init__(self, output_dir="alignment_analysis", min_frequency=0.01, top_n=10):
        """
        Initialize AlignmentAnalyzer
        
        Args:
            output_dir (str): Directory to store analysis results
            min_frequency (float): Minimum frequency threshold (0-1)
            top_n (int): Number of top sequences to keep per sample
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_frequency = min_frequency
        self.top_n = top_n
        
    def _extract_sequence_info(self, header, sequence):
        """Extract metadata from sequence header"""
        # Parse Illumina header format
        # @INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y READ:FILTERED:FLAG:INDEX
        parts = header[1:].split()  # Remove @ and split
        id_parts = parts[0].split(':')
        
        return {
            'id': parts[0],
            'instrument': id_parts[0],
            'run': id_parts[1],
            'flowcell': id_parts[2],
            'lane': id_parts[3],
            'sequence': sequence,
            'length': len(sequence),
            'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence)
        }
    
    def analyze_file(self, input_file):
        """
        Analyze sequences in a FASTQ file
        
        Args:
            input_file (str): Path to input FASTQ file
            
        Returns:
            tuple: (path to FASTA file, path to analysis report)
        """
        input_path = Path(input_file)
        fasta_path = self.output_dir / f"{input_path.stem}_frequent.fasta"
        report_path = self.output_dir / f"{input_path.stem}_analysis.json"
        stats_path = self.output_dir / f"{input_path.stem}_stats.txt"
        
        try:
            logger.info(f"Analyzing file: {input_path.name}")
            
            # Count sequence frequencies
            sequence_counter = Counter()
            sequence_info = {}
            total_sequences = 0
            
            logger.info("Counting sequence frequencies...")
            with open(input_file) as f:
                while True:
                    header = f.readline().strip()
                    if not header:
                        break
                    sequence = f.readline().strip()
                    f.readline()  # plus line
                    f.readline()  # quality
                    
                    total_sequences += 1
                    sequence_counter[sequence] += 1
                    
                    # Store info for first occurrence of each sequence
                    if sequence not in sequence_info:
                        sequence_info[sequence] = self._extract_sequence_info(header, sequence)
            
            # Find frequent sequences
            frequent_sequences = []
            for sequence, count in sequence_counter.most_common(self.top_n):
                frequency = count / total_sequences
                if frequency >= self.min_frequency:
                    info = sequence_info[sequence]
                    info['count'] = count
                    info['frequency'] = frequency
                    frequent_sequences.append(info)
            
            # Write frequent sequences to FASTA
            logger.info("Writing frequent sequences to FASTA...")
            with open(fasta_path, 'w') as f:
                for i, info in enumerate(frequent_sequences, 1):
                    seq_id = f"freq_{i}_count_{info['count']}"
                    f.write(f">{seq_id}\n{info['sequence']}\n")
            
            # Save detailed analysis
            analysis_report = {
                'input_file': input_path.name,
                'total_sequences': total_sequences,
                'unique_sequences': len(sequence_counter),
                'frequent_sequences': frequent_sequences,
                'analysis_parameters': {
                    'min_frequency': self.min_frequency,
                    'top_n': self.top_n
                }
            }
            
            with open(report_path, 'w') as f:
                json.dump(analysis_report, f, indent=2)
            
            # Save readable statistics
            with open(stats_path, 'w') as f:
                f.write("=== Sequence Analysis Statistics ===\n\n")
                f.write(f"Input file: {input_path.name}\n\n")
                f.write("Summary:\n")
                f.write(f"Total sequences analyzed: {total_sequences}\n")
                f.write(f"Unique sequences found: {len(sequence_counter)}\n")
                f.write(f"Sequences above {self.min_frequency*100}% frequency: {len(frequent_sequences)}\n\n")
                
                f.write("Most frequent sequences:\n")
                for i, info in enumerate(frequent_sequences, 1):
                    f.write(f"\n{i}. Sequence Info:\n")
                    f.write(f"   Count: {info['count']}\n")
                    f.write(f"   Frequency: {info['frequency']*100:.2f}%\n")
                    f.write(f"   Length: {info['length']}\n")
                    f.write(f"   GC Content: {info['gc_content']*100:.2f}%\n")
                    f.write(f"   Sequence: {info['sequence'][:50]}...")
            
            logger.info(f"\nAnalysis completed:")
            logger.info(f"- Total sequences: {total_sequences}")
            logger.info(f"- Unique sequences: {len(sequence_counter)}")
            logger.info(f"- Frequent sequences: {len(frequent_sequences)}")
            logger.info(f"- FASTA file: {fasta_path}")
            logger.info(f"- Analysis report: {report_path}")
            logger.info(f"- Statistics: {stats_path}")
            
            return str(fasta_path), str(report_path)
            
        except Exception as e:
            logger.error(f"Error during analysis: {str(e)}")
            raise
