import pysam
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
import logging
import argparse
import json

class CoverageAnalyzer:
    def __init__(self, bam_file, bed_file=None, reference=None):
        """Initialize the coverage analyzer with input files."""
        self.bam_file = bam_file
        self.bed_file = bed_file
        self.reference = reference
        self.setup_logging()
        
    def setup_logging(self):
        """Configure logging for the analysis."""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
    def validate_input(self):
        """Validate input BAM file and its index."""
        try:
            bam = pysam.AlignmentFile(self.bam_file, "rb")
            if not bam.has_index():
                self.logger.warning("BAM index not found, creating index...")
                pysam.index(self.bam_file)
            return True
        except Exception as e:
            self.logger.error(f"Error validating input: {str(e)}")
            return False
    
    def read_bed_file(self,bed_file):
        regions = []
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#'):  # Skip comments
                    continue
                parts = line.strip().split()
                regions.append({
                    'chrom': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2])
                })
        return regions

        
    def calculate_coverage_for_region(self, region):
        """Calculate coverage for a specific genomic region."""
        coverage_array = []
        try:
            with pysam.AlignmentFile(self.bam_file, "rb") as bam:
                # Calculate coverage for the region
                for pileup_column in bam.pileup(
                    region['chrom'],
                    region['start'],
                    region['end'],
                    truncate=True
                ):
                    #coverage_value = pileup_column.n if pileup_column.n is not None else 0
                    print(pileup_column.n)
                    coverage_array.append(pileup_column.n)

                    
            return {
                'region': region,
                'coverage': np.array(coverage_array),
                'mean_coverage': np.mean(coverage_array),
                'median_coverage': np.median(coverage_array),
                'std_coverage': np.std(coverage_array)
            }
        except Exception as e:
            self.logger.error(f"Error processing region {region}: {str(e)}")
            return None
    
    def analyze_gc_bias(self, coverage_results):
        """Analyze GC content bias in coverage."""
        gc_bias_data = defaultdict(list)
        if self.reference:
            with pysam.FastaFile(self.reference) as fasta:
                for result in coverage_results:
                    region = result['region']
                    sequence = fasta.fetch(
                        region['chrom'],
                        region['start'],
                        region['end']
                    )
                    gc_content = (
                        sequence.count('G') + sequence.count('C')
                    ) / len(sequence)
                    gc_bias_data['gc_content'].append(gc_content)
                    gc_bias_data['mean_coverage'].append(result['mean_coverage'])
        gc_bias= pd.DataFrame(gc_bias_data)
        gc_bias.to_csv("gc_bias.csv")
        return pd.DataFrame(gc_bias_data)            
    
    def process_regions_parallel(self, regions, num_processes=4):
        """Process multiple regions in parallel."""
        results = []
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            future_to_region = {
                executor.submit(self.calculate_coverage_for_region, region): region 
                for region in regions
            }
            for future in future_to_region:
                result = future.result()
                if result:
                    results.append(result)
        return results
        

    
    def generate_coverage_report(self, coverage_results):
    #"""Generate comprehensive coverage report as a DataFrame."""
    # Calculate global statistics
        global_stats = {
            'mean_coverage': np.mean([r['mean_coverage'] for r in coverage_results]),
            'median_coverage': np.median([r['median_coverage'] for r in coverage_results]),
            'std_coverage': np.mean([r['std_coverage'] for r in coverage_results])
        }

        # Create DataFrame for global stats
        global_stats_df = pd.DataFrame(global_stats, index=[0])
        global_stats_df.to_csv("global_stats.csv")

        # Create DataFrame for regions
        regions_data = [{
            'chrom': r['region']['chrom'],
            'start': r['region']['start'],
            'end': r['region']['end'],
            'mean_coverage': r['mean_coverage'],
            'median_coverage': r['median_coverage'],
            'std_coverage': r['std_coverage']
        } for r in coverage_results]

        regions_df = pd.DataFrame(regions_data)
        regions_df.to_csv("region.csv")
        

        # Calculate coverage distribution percentiles
        all_coverage = np.concatenate([r['coverage'] for r in coverage_results])
        regions_low = regions_df[regions_df["mean_coverage"] < np.percentile(all_coverage, 30)] 
        regions_low.to_csv("regions_low.csv")
        regions_high = regions_df[regions_df["mean_coverage"] > np.percentile(all_coverage, 70)] 
        regions_high.to_csv("regions_high.csv")
        coverage_table = pd.DataFrame({
            'position': range(len(all_coverage)),
            'depth': all_coverage
        })
        coverage_table.to_csv("coverage_depth.csv")

        percentiles = {
            'p10': np.percentile(all_coverage, 10),
            'p25': np.percentile(all_coverage, 25),
            'p50': np.percentile(all_coverage, 50),
            'p75': np.percentile(all_coverage, 75),
            'p90': np.percentile(all_coverage, 90)
        }
        
        # Create DataFrame for percentiles
        percentiles_df = pd.DataFrame(percentiles, index=[0])
        percentiles_df.to_csv("percentiles.csv")

        # Concatenate all DataFrames into one final report
        report_df = pd.concat([global_stats_df, percentiles_df, regions_df], axis=0, ignore_index=True)

        return report_df

        
    def plot_coverage_distribution(self, coverage_results, output_file):
        
        all_coverage = np.concatenate([r['coverage'] for r in coverage_results])
        
        plt.figure(figsize=(10, 6))
        plt.hist(all_coverage, bins=1000, alpha=0.75)
        plt.xlim(0, 20)
        plt.xlabel('Coverage Depth')
        plt.ylabel('Frequency')
        plt.title('Coverage Distribution')
        plt.axvline(np.mean(all_coverage), color='r', linestyle='dashed',
                   label=f'Mean Coverage: {np.mean(all_coverage):.2f}')
        plt.legend()
        plt.savefig(output_file)
        plt.close()
    

    
def main():
    parser = argparse.ArgumentParser(description='NGS Coverage Analysis Tool')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--bed', help='BED file with regions of interest')
    parser.add_argument('--reference', help='Reference genome FASTA')
    parser.add_argument('--output', default='coverage_report',
                       help='Output prefix for results')
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = CoverageAnalyzer(args.bam, args.bed, args.reference)
    
    # Validate input
    if not analyzer.validate_input():
        return
    
    # Define regions (from BED file or whole genome)
    #regions = [{'chrom': 'chr1', 'start': 0, 'end': 1000000}]  # Example region
    if args.bed:
        print(args.bed)
        regions = analyzer.read_bed_file(args.bed)
    else:
        regions = [{'chrom': 'chr1', 'start': 0, 'end': 1000000}]  # Example region
    # Calculate coverage
    coverage_results = analyzer.process_regions_parallel(regions)
    
    # Generate report
    report_df = analyzer.generate_coverage_report(coverage_results)
    
    
    #generate gc_bias table
    analyzer.analyze_gc_bias(coverage_results)
    # Plot distribution
    analyzer.plot_coverage_distribution(
        coverage_results,
        f"{args.output}_distribution.png"
    )
  
if __name__ == "__main__":
    main()
