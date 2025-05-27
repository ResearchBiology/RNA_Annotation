import logging
from pathlib import Path
import subprocess
import json
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time
from tqdm import tqdm

logger = logging.getLogger(__name__)

class BlastSearcher:
    """Class for performing BLAST searches on sequences"""
    
    def __init__(self, output_dir="blast_results", blast_db=None, email=None):
        """
        Initialize BlastSearcher
        
        Args:
            output_dir (str): Directory to store BLAST results
            blast_db (str): Path to local BLAST database (optional)
            email (str): Email for NCBI BLAST (required for online BLAST)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.blast_db = blast_db
        self.email = email
    
    def _run_local_blast(self, input_file, output_file):
        """Run BLAST search using local database"""
        cmd = [
            'blastn',
            '-query', str(input_file),
            '-db', self.blast_db,
            '-out', str(output_file),
            '-outfmt', '15',  # JSON format
            '-max_target_seqs', '5',
            '-evalue', '1e-10'
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST error: {e.stderr}")
            return False
    
    def _run_online_blast(self, sequence_record, output_file):
        """Run BLAST search using NCBI web service"""
        try:
            # Use XML format as it's more stable than JSON for NCBI BLAST
            result_handle = NCBIWWW.qblast(
                "blastn",
                "nt",
                sequence_record.format("fasta"),
                hitlist_size=5,
                expect=1e-10
            )
            
            with open(output_file, 'w') as f:
                f.write(result_handle.read())
            
            result_handle.close()
            return True
            
        except Exception as e:
            logger.error(f"NCBI BLAST error: {str(e)}")
            return False
    
    def _parse_blast_results(self, blast_file):
        """Parse BLAST results from XML format"""
        hits = []
        with open(blast_file) as f:
            blast_records = NCBIXML.parse(f)
            for blast_record in blast_records:
                for alignment in blast_record.alignments[:5]:  # Top 5 hits
                    hsp = alignment.hsps[0]  # Best HSP
                    
                    # Extract accession and species from title
                    title_parts = alignment.title.split('|')
                    accession = title_parts[1] if len(title_parts) > 1 else ''
                    description = title_parts[-1].strip()
                    
                    # Try to extract species name
                    species = ''
                    if '[' in description and ']' in description:
                        species = description[description.rfind('[')+1:description.rfind(']')]
                    
                    hits.append({
                        'title': description,
                        'accession': accession,
                        'scientific_name': species,
                        'percent_identity': (hsp.identities / hsp.align_length) * 100,
                        'alignment_length': hsp.align_length,
                        'e_value': hsp.expect,
                        'bit_score': hsp.bits
                    })
        
        return hits
    
    def search(self, input_file):
        """
        Perform BLAST search on sequences in FASTA file
        
        Args:
            input_file (str): Path to input FASTA file
            
        Returns:
            str: Path to results file
        """
        input_path = Path(input_file)
        results_file = self.output_dir / f"{input_path.stem}_blast_results.json"
        summary_file = self.output_dir / f"{input_path.stem}_blast_summary.txt"
        
        try:
            logger.info(f"Processing file: {input_path.name}")
            all_results = []
            
            # Read sequences
            sequences = list(SeqIO.parse(input_file, "fasta"))
            logger.info(f"Found {len(sequences)} sequences to search")
            
            for seq_record in tqdm(sequences, desc="Running BLAST searches"):
                seq_file = self.output_dir / f"temp_{seq_record.id}.fasta"
                blast_file = self.output_dir / f"temp_{seq_record.id}_blast.json"
                
                # Write individual sequence to temp file
                SeqIO.write(seq_record, seq_file, "fasta")
                
                # Run BLAST
                success = False
                if self.blast_db:
                    logger.info(f"Running local BLAST for {seq_record.id}")
                    success = self._run_local_blast(seq_file, blast_file)
                elif self.email:
                    logger.info(f"Running online BLAST for {seq_record.id}")
                    success = self._run_online_blast(seq_record, blast_file)
                else:
                    raise ValueError("Neither local BLAST database nor email for NCBI BLAST provided")
                
                if success:
                    # Parse results
                    hits = self._parse_blast_results(blast_file)
                    all_results.append({
                        'query_id': seq_record.id,
                        'query_length': len(seq_record.seq),
                        'hits': hits
                    })
                
                # Clean up temp files
                seq_file.unlink(missing_ok=True)
                blast_file.unlink(missing_ok=True)
                
                # Wait between online BLAST requests
                if self.email and not self.blast_db:
                    time.sleep(60)  # NCBI requires delay between requests
            
            # Save all results
            with open(results_file, 'w') as f:
                json.dump(all_results, f, indent=2)
            
            # Create human-readable summary
            with open(summary_file, 'w') as f:
                f.write("=== BLAST Search Results ===\n\n")
                f.write(f"Input file: {input_path.name}\n")
                f.write(f"Total sequences searched: {len(sequences)}\n\n")
                
                for result in all_results:
                    f.write(f"\nQuery: {result['query_id']} (Length: {result['query_length']})\n")
                    if result['hits']:
                        f.write("Top hits:\n")
                        for i, hit in enumerate(result['hits'], 1):
                            f.write(f"\n{i}. {hit['title']}\n")
                            f.write(f"   Accession: {hit['accession']}\n")
                            f.write(f"   Species: {hit['scientific_name']}\n")
                            f.write(f"   Identity: {hit['percent_identity']:.2f}%\n")
                            f.write(f"   E-value: {hit['e_value']:.2e}\n")
                            f.write(f"   Bit score: {hit['bit_score']:.1f}\n")
                    else:
                        f.write("No significant hits found\n")
            
            logger.info(f"\nBLAST search completed:")
            logger.info(f"- Results saved to: {results_file}")
            logger.info(f"- Summary saved to: {summary_file}")
            
            return str(results_file)
            
        except Exception as e:
            logger.error(f"Error during BLAST search: {str(e)}")
            raise
