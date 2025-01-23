import logging
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

logger = logging.getLogger(__name__)

class QualityFilter:
    """Filter sequences based on quality scores using FastQC-style analysis"""
    
    def __init__(self, output_dir="results/quality", window_size=15, min_quality=30):
        """
        Initialize QualityFilter
        
        Args:
            output_dir (str): Directory to store filtered sequences
            window_size (int): Size of sliding window for quality calculation
            min_quality (int): Minimum quality score (Q30 by default)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.window_size = window_size
        self.min_quality = min_quality
        self.plots_dir = self.output_dir / "plots"
        self.plots_dir.mkdir(exist_ok=True)
    
    def _calculate_quality_scores(self, quality_string):
        """Convert ASCII quality string to numeric scores"""
        return [ord(char) - 33 for char in quality_string]  # Phred+33 encoding
    
    def _check_window_quality(self, quality_scores):
        """Check if any window falls below minimum quality"""
        for i in range(len(quality_scores) - self.window_size + 1):
            window = quality_scores[i:i + self.window_size]
            if np.mean(window) < self.min_quality:
                return False
        return True
    
    def _generate_quality_plot(self, position_qualities, output_path):
        """Generate quality score distribution plot"""
        plt.figure(figsize=(12, 6))
        
        # Calculate statistics for each position
        positions = sorted(position_qualities.keys())
        medians = [np.median(position_qualities[p]) for p in positions]
        q1 = [np.percentile(position_qualities[p], 25) for p in positions]
        q3 = [np.percentile(position_qualities[p], 75) for p in positions]
        
        # Plot quality distribution
        plt.fill_between(positions, q1, q3, alpha=0.2, color='blue')
        plt.plot(positions, medians, color='blue', label='Median Score')
        
        # Add threshold line
        plt.axhline(y=self.min_quality, color='red', linestyle='--', 
                   label=f'Q{self.min_quality} Threshold')
        
        plt.title('Quality Scores Distribution Across All Positions')
        plt.xlabel('Position in Read (bp)')
        plt.ylabel('Quality Score')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Save plot
        plt.savefig(output_path)
        plt.close()
    
    def filter_sequences(self, input_file):
        """
        Filter sequences based on quality scores
        
        Args:
            input_file (str): Path to input file
            
        Returns:
            str: Path to filtered file
        """
        try:
            input_path = Path(input_file)
            output_path = self.output_dir / f"quality_filtered_{input_path.name}"
            stats_path = self.output_dir / f"{input_path.stem}_quality_stats.txt"
            plot_path = self.plots_dir / f"{input_path.stem}_quality_plot.png"
            
            # Track statistics
            total_sequences = 0
            passed_sequences = 0
            position_qualities = defaultdict(list)
            
            logger.info(f"Processing file: {input_path.name}")
            
            # Process file
            with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
                while True:
                    # Read sequence entry (2 lines)
                    info_line = infile.readline().strip()
                    seq_line = infile.readline().strip()
                    
                    # Check for end of file
                    if not info_line or not seq_line:
                        break
                    
                    total_sequences += 1
                    
                    # Convert quality scores
                    quality_scores = self._calculate_quality_scores(seq_line)
                    
                    # Track qualities for plotting
                    for pos, score in enumerate(quality_scores):
                        position_qualities[pos].append(score)
                    
                    # Check quality using sliding window
                    if self._check_window_quality(quality_scores):
                        outfile.write(f"{info_line}\n{seq_line}\n")
                        passed_sequences += 1
                    
                    # Log progress
                    if total_sequences % 100000 == 0:
                        logger.info(f"Processed {total_sequences} sequences")
            
            # Generate quality plot
            self._generate_quality_plot(position_qualities, plot_path)
            
            # Calculate statistics
            failed_sequences = total_sequences - passed_sequences
            
            # Save statistics
            with open(stats_path, 'w') as f:
                f.write("=== Quality Filtering Statistics ===\n\n")
                f.write(f"Input file: {input_path.name}\n")
                f.write(f"Quality threshold: Q{self.min_quality}\n")
                f.write(f"Window size: {self.window_size}\n\n")
                f.write("Summary:\n")
                f.write(f"Total sequences: {total_sequences}\n")
                f.write(f"Sequences passed: {passed_sequences} ({passed_sequences/total_sequences*100:.2f}%)\n")
                f.write(f"Sequences failed: {failed_sequences} ({failed_sequences/total_sequences*100:.2f}%)\n")
                
                # Add quality score distribution
                f.write("\nQuality Score Distribution:\n")
                for pos in sorted(position_qualities.keys()):
                    scores = position_qualities[pos]
                    f.write(f"Position {pos}: Mean={np.mean(scores):.2f}, "
                           f"Median={np.median(scores):.2f}, "
                           f"Min={min(scores)}, Max={max(scores)}\n")
            
            logger.info(f"\nQuality filtering completed:")
            logger.info(f"- Total sequences: {total_sequences}")
            logger.info(f"- Sequences passed: {passed_sequences}")
            logger.info(f"- Sequences failed: {failed_sequences}")
            logger.info(f"- Results saved to: {output_path}")
            logger.info(f"- Statistics saved to: {stats_path}")
            logger.info(f"- Quality plot saved to: {plot_path}")
            
            return str(output_path)
            
        except Exception as e:
            logger.error(f"Error processing file {input_file}: {str(e)}")
            raise 