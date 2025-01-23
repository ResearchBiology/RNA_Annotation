import logging
from pathlib import Path
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

logger = logging.getLogger(__name__)

class BlastPDFGenerator:
    """Class for generating PDF reports from BLAST results"""
    
    def __init__(self, blast_dir="results/blast", output_dir="results/reports"):
        """Initialize PDF Generator"""
        self.blast_dir = Path(blast_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.styles = getSampleStyleSheet()
        
        # Create custom styles
        self.styles.add(ParagraphStyle(
            name='CustomTitle',
            parent=self.styles['Heading1'],
            fontSize=14,
            spaceAfter=30
        ))
        self.styles.add(ParagraphStyle(
            name='CustomBody',
            parent=self.styles['Normal'],
            fontSize=10,
            spaceAfter=12
        ))
    
    def _parse_blast_file(self, file_path):
        """Parse a BLAST result file and return structured data"""
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Extract sample name from first line
        sample_name = lines[0].strip().replace("BLAST Results for ", "")
        
        # Find query length
        query_length = next(line for line in lines if "Query Length:" in line)
        
        # Parse hits
        hits = []
        current_hit = []
        in_hit = False
        
        for line in lines:
            if line.startswith("Hit #"):
                if current_hit:
                    hits.append(current_hit)
                current_hit = [line.strip()]
                in_hit = True
            elif in_hit and line.strip() and not line.startswith("-"):
                current_hit.append(line.strip())
        
        if current_hit:
            hits.append(current_hit)
            
        return {
            'sample_name': sample_name,
            'query_length': query_length.strip(),
            'hits': hits
        }
    
    def _create_hit_table(self, hit_data):
        """Create a formatted table for a BLAST hit"""
        data = []
        for line in hit_data:
            data.append([Paragraph(line, self.styles['CustomBody'])])
        
        table = Table(data, colWidths=[6.5*inch])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, 0), colors.lightgrey),
            ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 12),
            ('TOPPADDING', (0, 0), (-1, -1), 12),
            ('BOX', (0, 0), (-1, -1), 1, colors.black),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey)
        ]))
        return table
    
    def generate_report(self):
        """Generate a PDF report containing all BLAST results"""
        try:
            # Get all BLAST result files
            blast_files = sorted(self.blast_dir.glob("*_blast.txt"))
            if not blast_files:
                logger.warning("No BLAST result files found")
                return None
            
            # Create PDF document
            output_file = self.output_dir / "blast_results_report.pdf"
            doc = SimpleDocTemplate(
                str(output_file),
                pagesize=letter,
                rightMargin=72,
                leftMargin=72,
                topMargin=72,
                bottomMargin=72
            )
            
            # Container for PDF elements
            elements = []
            
            # Process each BLAST file
            for blast_file in blast_files:
                # Parse BLAST results
                result_data = self._parse_blast_file(blast_file)
                
                # Add sample title
                elements.append(Paragraph(
                    f"Sample: {result_data['sample_name']}", 
                    self.styles['CustomTitle']
                ))
                
                # Add query length
                elements.append(Paragraph(
                    result_data['query_length'],
                    self.styles['CustomBody']
                ))
                elements.append(Spacer(1, 0.2*inch))
                
                # Add each hit
                for hit in result_data['hits']:
                    elements.append(self._create_hit_table(hit))
                    elements.append(Spacer(1, 0.2*inch))
                
                # Add page break
                elements.append(Paragraph("<br clear=all style='page-break-before:always'/>", 
                                       self.styles['CustomBody']))
            
            # Build PDF
            doc.build(elements)
            logger.info(f"PDF report generated: {output_file}")
            return str(output_file)
            
        except Exception as e:
            logger.error(f"Error generating PDF report: {str(e)}")
            return None 