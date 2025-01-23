import logging
from pathlib import Path
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time

logger = logging.getLogger(__name__)

def process_single_file(file_path, output_dir):
    """Process a single file for BLAST search"""
    try:
        logger.info(f"Processing BLAST search for {file_path.name}")
        
        # Get sample name
        sample_name = file_path.stem.replace("_sequence", "")
        
        # Read sequence
        with open(file_path, 'r') as f:
            lines = f.readlines()
            sequence = None
            for i, line in enumerate(lines):
                if line.startswith("Sequence:"):
                    sequence = lines[i+1].strip()
                    break
        
        if not sequence:
            logger.warning(f"No sequence found in {file_path.name}")
            return None
        
        # Add delay before request
        time.sleep(5)
        
        # Perform simple BLAST search
        logger.info(f"Starting BLAST search for {sample_name}")
        result_handle = NCBIWWW.qblast(
            "blastn",
            "nt",
            sequence
        )
        
        # Parse results
        blast_records = NCBIXML.parse(result_handle)
        record = next(blast_records)
        
        # Write results
        output_file = Path(output_dir) / f"{sample_name}_blast.txt"
        with open(output_file, 'w') as f:
            f.write(f"BLAST Results for {sample_name}\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Search Type: BLASTN\n")
            f.write(f"Query Length: {record.query_length} bp\n\n")
            f.write("Top Hits:\n")
            f.write("-" * 80 + "\n")
            
            for i, alignment in enumerate(record.alignments, 1):
                hsp = alignment.hsps[0]
                f.write(f"\nHit #{i}:\n")
                f.write(f"Title: {alignment.title}\n")
                f.write(f"Length: {alignment.length} bp\n")
                f.write(f"E-value: {hsp.expect:.2e}\n")
                f.write(f"Identity: {hsp.identities}/{hsp.align_length} ({(hsp.identities/hsp.align_length)*100:.1f}%)\n")
                f.write(f"Gaps: {hsp.gaps}/{hsp.align_length} ({(hsp.gaps/hsp.align_length)*100:.1f}%)\n")
                f.write(f"Query: {hsp.query[0:50]}...\n")
                f.write(f"Match: {hsp.match[0:50]}...\n")
                f.write(f"Sbjct: {hsp.sbjct[0:50]}...\n")
                f.write("-" * 40 + "\n")
        
        logger.info(f"Completed BLAST search for {sample_name}")
        return str(output_file)
        
    except Exception as e:
        logger.error(f"Error processing BLAST search for {file_path.name}: {str(e)}")
        return None

class BlastSearcher:
    """Class for performing BLAST searches on sequences"""
    
    def __init__(self, input_dir="results/alignment", output_dir="results/blast"):
        """Initialize BlastSearcher"""
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def process_files(self):
        """Process all alignment files and perform BLAST searches sequentially"""
        # Get list of files
        input_files = sorted(self.input_dir.glob("*_sequence.txt"))
        processed_files = []
        
        # Process files one at a time
        for input_file in input_files:
            result = process_single_file(input_file, self.output_dir)
            if result:
                processed_files.append(result)
                # Add delay between files
                time.sleep(45)
        
        logger.info(f"\nBLAST Processing completed:")
        logger.info(f"- Total files processed: {len(processed_files)}")
        logger.info(f"- Results saved in: {self.output_dir}")
        
        return processed_files 