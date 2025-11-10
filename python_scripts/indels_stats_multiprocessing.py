import os
import glob
import re
import sys
import pysam
import pandas as pd
# multiprocessing packages
from multiprocessing import Pool

# module load statements in unity
'''
module load uri/main
module load Pysam/0.17.0-GCC-11.2.0

'''

# process one sample VCF and return variant data/stats
def ProcessSample(vcf_path):
    sample = os.path.basename(vcf_path).split('_')[0] # extract the sample from the file name
    summary = [] #output 
        
    with pysam.VariantFile(vcf_path) as vcf:
        for entry in vcf.fetch():
            chrom = entry.chrom     # extract the chromosome
            pos = entry.pos         # extract the variant position
            var_id = entry.id or f"{chrom}_{pos}"  # extract variant ID or create one with the chromosome and possition if missing

            info = entry.info       # extract the info field from the annotations for each variant
            
            # calculate the length of the sv
            end = info.get("END", pos + len(entry.ref)-1) # extracts the end of the variant or estimates based on reference allele len - 1 if missing
            svlen = abs(info.get("SVLEN", end-pos)) # extracts the length or calculates using the end and pos if missing
            
            svtype = info.get("SVTYPE", "NA")   # extracts the type (insertion or deletion) or places NA if it doesn't exist
            
            # extracting the genes associated with annotation(s) attached to each variant
            annotations = entry.info.get("ANN", "")
            if annotations:
                gene_ids = []
                for ann in annotations.split(","):  # if multiple annotations, they are split with a comma
                    fields = ann.split("|")         # each field in the annotations is split by a |
                    if len(fields) > 4 and fields[4]:
                        gene_ids.append(fields[4])  # add each gene_id
                gene_ids_unique = list(set(gene_ids))  # remove any duplicate gene ids
            else:
                gene_ids_unique = []
            
            # add all info from each variant to the output
            summary.append({
                "Sample": sample,
                "Variant_ID": var_id,
                "Chromosome": chrom,
                "Start": pos,
                "End": end,
                "Size": svlen,
                "Variant_Type": svtype,
                "Gene_IDs": gene_ids_unique
            })
    return summary
                
if __name__ == "__main__":
    base_dir = sys.argv[1]      # in this case thats "/work/pi_nhowlett_uri_edu/jessie/New-All-20-Bam/""
    output_path = sys.argv[2]   # path to where the output tsv should go
    num_samples = int(sys.argv[3]) # input number of samples (n=20)
    
    # extract all *ann.vcf files from the base dir
    vcf_files = glob.glob(os.path.join(base_dir, "[0-9]*-bc*/gatk/*.ann.vcf"))
   
    print(f"Found {len(vcf_files)} files:")
    for f in vcf_files:
        print(f)

    print(f"Processing {len(vcf_files)} files using {num_samples} processes")
    
    with Pool(processes=num_samples) as pool:
        results = pool.map(ProcessSample, vcf_files)
    
    # flatten the list of lists that we get from the multiprocessing
    all_records = [variant for sample_results in results for variant in sample_results]
    
    # write the output to a tsv file
    df = pd.DataFrame(all_records)
    df.to_csv(output_path, sep="\t", index=False) # writes to the tsv path specified
    print("\n Saved parsed variants to {output_path}")
