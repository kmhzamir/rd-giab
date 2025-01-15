process MERGE_VCF_COLUMNS {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/kmhzamir/dragena:v1.4"
    
    input:
    tuple val(meta), path(merged_tab)  // Concatenated tab file from previous step
    tuple val(meta), path(vcf), path(tbi)                           // Raw VCF file

    output:
    tuple val(meta), path("${meta.id}_mt_final.tab.gz")   , emit: cleaned_output

    script:
    """
    python3 <<CODE
import gzip
import pandas as pd

# Define the file paths
vcf_file_path = '${vcf}'
annotated_file_path = '${merged_tab}'
output_file_path = '${meta.id}_merged_with_vcf_data.tab.gz'

# Load the VCF file and dynamically get the 10th column name
vcf_columns = []
vcf_data = []
with gzip.open(vcf_file_path, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            vcf_columns = line.strip().split('\\t')
            sample_name = vcf_columns[9]  # Get the name of the 10th column dynamically
            break

# Now read the VCF data lines and extract relevant columns
with gzip.open(vcf_file_path, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue  # Skip header lines
        fields = line.strip().split('\\t')
        vcf_data.append({
            '#CHROM': fields[0],
            'POS': fields[1],
            'QUAL': fields[5],
            'INFO': fields[7],
            'FORMAT': fields[8],
            sample_name: fields[9]
        })

# Convert to DataFrame
vcf_df = pd.DataFrame(vcf_data, columns=['#CHROM', 'POS', 'QUAL', 'INFO', 'FORMAT', sample_name])
# Create Location column in format '#CHROM:POS'
vcf_df['Location'] = vcf_df['#CHROM'] + ':' + vcf_df['POS']

# Find and load the annotated file header row
with gzip.open(annotated_file_path, 'rt') as f:
    header_row = 0
    for line in f:
        if line.startswith('#') and not line.startswith('##'):
            header = line.strip().split('\\t')  # Keep this header line
            break
        header_row += 1

# Now load the data from the annotated file, skipping all lines up to the found header
annotated_df = pd.read_csv(
    annotated_file_path,
    sep='\\t',
    compression='gzip',
    header=header_row
)
annotated_df.columns = header  # Set the header correctly

# Check if 'Location' exists in annotated_df for merging
if 'Location' not in annotated_df.columns:
    print("The 'Location' column is not found in annotated_df. Please check the column names.")
    exit(1)

# Merge data based on the Location column
merged_df = pd.merge(annotated_df, vcf_df[['Location', 'QUAL', 'INFO', 'FORMAT', sample_name]], 
                     on='Location', how='left')

# Save the result to a new file
merged_df.to_csv(output_file_path, sep='\\t', index=False, compression='gzip')
print("Merged file saved to:", output_file_path)
CODE

    # Run zcat, grep, awk, and bgzip to clean up duplicates in the merged file
    zcat ${meta.id}_merged_with_vcf_data.tab.gz | awk '!/^#/ {print}' | bgzip > ${meta.id}_mt_final.tab.gz
    """
}
