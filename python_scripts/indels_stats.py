import os
import statistics
import sys
from collections import defaultdict
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# multiprocessing packages
from multiprocessing import Pool

def GetGeneAccession(variant, info_elements):
    accession = 'NA'
    # get gene(s) associated with variant
    for e in info_elements:
        if e.split('=')[0] == 'ANN':
            annotations = e.split(',')
            annotations[0] = annotations[0][4:]
            for a in annotations:
                # get variant
                ann_v = a.split('|')[0]
                if ann_v == variant:
                    if a.split('|')[5] in ['transcript', 'protein_coding']:
                        accession = a.split('|')[6]
    return accession

# process one sample VCF and return variant data/stats
def ProcessSample(args):
    vcf_file, sample, output_dir = args
    
    # store the size of indels, number of multiple variants and number of variants
    insertions = []
    deletions = []
    mul_variants = 0
    num_variants = 0
    summary_records = []
    
    with open(vcf_file, 'r') as f:
        for count, line in enumerate(f, 1):
            # ignore header (lines that start with #)
            if line.rstrip()[0] != '#':
                # get reference and variant
                ref = line.rstrip().split('\t')[3]
                variants = line.rstrip().split('\t')[4]
                info_elements = line.rstrip().split('\t')[7].split(';')
                for v in variants:
                    v_size = int()
                    v_type = ''
                    # identify if variant is an insertion or a deletion
                    if len(ref) > len(v):
                        v_size = len(ref) - len(v)
                        v_type = 'deletion'
                    if len(v) > len(ref):
                        v_size = len(v) - len(ref)
                        v_type = 'insertion'
                        
                    if v_size >= variant_size_threshold:
                        num_variants += 1
                        accession = GetGeneAccession(v, info_elements)
                        if 'unassigned_transcript_' in accession:
                            accession = 'NA'
                        # gather information about each indel and add to the summary records
                        if v_type == 'insertion':
                            record = {
                                'Sample': sample,
                                'Variant': v, 
                                'Reference': ref,
                                'Size': v_size,
                                'Accession': accession,
                                'Type': 'insertion'
                            }
                            summary_records.append(record)
                            insertions.append(v_size)
                        elif v_type == 'deletion':
                            record = {
                                'Sample': sample,
                                'Variant': v, 
                                'Reference': ref,
                                'Size': v_size,
                                'Accession': accession,
                                'Type': 'deletion'
                            }
                            summary_records.append(record)  
                            deletions.append(v_size)
                if len(variants) > 1:
                    mul_variants += 1
    return summary_records


if __name__ == "__main__":
    vcf_dir = sys.argv[1]
    output_dir = sys.argv[2]
    num_samples = int(sys.argv[3])
    variant_size_threshold = int(sys.argv[4])

    # creates full paths for each vcf file in the specified directory
    vcf_files = [os.path.join(vcf_dir, f) for f in os.listdir(vcf_dir) if f.endswith('.vcf')]
    # creates the triple variable for each sample's vcf for input into the process sample function
    tasks = [(vcf, os.path.basename(vcf).replace('.vcf', ''), output_dir) for vcf in vcf_files]

    # use the multiprocessing library tools to run in parallel
    print(f"Runing {len(vcf_files)} samples using {num_samples} processes")
    with Pool(processes=num_samples) as pool:
        results = pool.map(ProcessSample, tasks)
    
    all_records = []
    for records in results:
        all_records.extend(records)
    
    # create combined summary
    summary_path = os.path.join(output_dir, 'all_samples_cnv_info.tsv')
    df = pd.DataFrame(all_records)
    df.to_csv(summary_path, sep='\t', index=False)