import pandas as pd
import os
import glob
import re

def combine_ploidy_results(input_base_dir, output_dir):
    """
    Concatenates GATK CNV contig ploidy results into a single TSV table.

    Args:
        input_base_dir (str): Path to the base directory containing SAMPLE_X subdirectories.
        output_dir (str): Path to the directory where the output TSV will be saved.
    """
    all_ploidy_data = []

    # Find all SAMPLE_X directories
    sample_dirs = glob.glob(os.path.join(input_base_dir, "SAMPLE_*"))
    # Sort sample directories to ensure consistent processing order
    sample_dirs.sort(key=lambda x: int(re.search(r'SAMPLE_(\d+)', x).group(1)))


    if not sample_dirs:
        print(f"No 'SAMPLE_*' directories found in {input_base_dir}")
        return

    print(f"Found {len(sample_dirs)} sample directories.")

    for sample_dir in sample_dirs:
        sample_id = os.path.basename(sample_dir)
        
        # Read sample name
        sample_name_file = os.path.join(sample_dir, "sample_name.txt")
        if not os.path.exists(sample_name_file):
            print(f"Warning: {sample_name_file} not found. Skipping {sample_id}.")
            continue
        with open(sample_name_file, 'r') as f:
            strain_name = f.readline().strip()

        # Read contig ploidy
        contig_ploidy_file = os.path.join(sample_dir, "contig_ploidy.tsv")
        if not os.path.exists(contig_ploidy_file):
            print(f"Warning: {contig_ploidy_file} not found. Skipping {sample_id}.")
            continue

        # Read TSV, skipping lines starting with '@'
        try:
            df_ploidy = pd.read_csv(contig_ploidy_file, sep='	', comment='@')
            # Select only CONTIG and PLOIDY columns
            df_ploidy = df_ploidy[['CONTIG', 'PLOIDY']]
            df_ploidy['STRAIN'] = strain_name
            all_ploidy_data.append(df_ploidy)
        except Exception as e:
            print(f"Error reading {contig_ploidy_file}: {e}. Skipping {sample_id}.")
            continue

    if not all_ploidy_data:
        print("No ploidy data could be processed.")
        return

    # Concatenate all dataframes
    combined_df = pd.concat(all_ploidy_data, ignore_index=True)

    # Pivot the table: rows=strains, columns=chromosome, values=ploidy
    # Use reset_index() to make 'STRAIN' a regular column for pivoting
    pivot_table = combined_df.pivot(index='STRAIN', columns='CONTIG', values='PLOIDY')

    # Reorder columns to ensure chrI, chrII, ... order
    # Extract chromosome numbers and sort them numerically
    def sort_chromosomes_key(col_name):
        if col_name.startswith('chr'):
            # Extract Roman numeral and convert to integer for proper sorting
            roman_numeral = col_name[3:]
            # Simple mapping for common Roman numerals up to XVI
            roman_map = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 
                         'VIII': 8, 'IX': 9, 'X': 10, 'XI': 11, 'XII': 12, 'XIII': 13, 
                         'XIV': 14, 'XV': 15, 'XVI': 16}
            return roman_map.get(roman_numeral, 999) # Use a high number for unrecognised
        return 999 # Non-'chr' columns go to the end

    # Sort columns by chromosome number using the custom key
    sorted_columns = sorted(pivot_table.columns, key=sort_chromosomes_key)
    pivot_table = pivot_table[sorted_columns]

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    output_file = os.path.join(output_dir, "ploidy_table.tsv")
    pivot_table.to_csv(output_file, sep='\t')
    print(f"Ploidy table saved to: {output_file}")

    # --- Generate Abnormal Ploidy Summary ---
    abnormal_ploidy_summary = []
    for index, row in pivot_table.iterrows():
        abnormal_chromosomes_list = []
        abnormal_ploidy_values = set()
        for col in pivot_table.columns:
            if row[col] != 2:
                abnormal_chromosomes_list.append(col)
                abnormal_ploidy_values.add(int(row[col]))
        
        if abnormal_chromosomes_list:
            # Sort chromosomes for consistent output
            sorted_abnormal_chromosomes = sorted(abnormal_chromosomes_list, key=sort_chromosomes_key)
            # Sort ploidy values for consistent output
            sorted_abnormal_ploidy_values = sorted(list(abnormal_ploidy_values))

            abnormal_ploidy_summary.append({
                'STRAIN': index,
                'AbnormalChromosomes': ",".join(sorted_abnormal_chromosomes),
                'AbnormalPloidyValues': ",".join(map(str, sorted_abnormal_ploidy_values))
            })
    
    if abnormal_ploidy_summary:
        df_abnormal_ploidy = pd.DataFrame(abnormal_ploidy_summary)
        abnormal_output_file = os.path.join(output_dir, "abnormal_ploidy_summary.tsv")
        df_abnormal_ploidy.to_csv(abnormal_output_file, sep='\t', index=False)
        print(f"Abnormal ploidy summary saved to: {abnormal_output_file}")
    else:
        print("No abnormal ploidy detected in any sample.")


if __name__ == "__main__":
    input_base_dir = "gatk_cnv/ploidy-calls"
    output_dir = "tables"
    combine_ploidy_results(input_base_dir, output_dir)