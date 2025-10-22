import gzip
from collections import defaultdict
import math
import os

# Function to clean the individual ID by removing column numbers and ':GT'
def clean_individual_id(individual_id):
    return individual_id.split(":")[0].split("]")[-1]  # Remove column numbers and ':GT'

# Function to find the midpoint between two positions and round accordingly
def midpoint(start, end, round_up=True):
    mid = (start + end) / 2
    return math.ceil(mid) if round_up else math.floor(mid)

# Function to count consecutive homozygous alleles ("1/1" or "1|1") by row with header handling
def count_consecutive_homozygous_with_header(input_file, counting_genotypes, missing_genotypes):
    result = []
    header = None  # Placeholder for the header

    with gzip.open(input_file, 'rt') as file:
        
        for idx, line in enumerate(file):
            line = line.strip().split("\t")
            
            
            if idx == 0: #Treat the first line as the header
                header = line
                result.append(header)  # Append the header to the result as is
                current_chromosome = None
            else: #content
                chromosome, position = line[0], line[1]  # Keep chromosome and position
                genotypes = line[2:]  # Get all genotype data
                if current_chromosome is None: #for start of the first chromosome
                  previous_row_counts = None
                  current_chromosome = chromosome
                elif current_chromosome != chromosome: #for change to the next chromosome
                  previous_row_counts = [ 0 for i in previous_row_counts]
                  current_chromosome = chromosome
                new_row = [chromosome, position]  # Start the result row with chromosome and position
                current_counts = []
                
                # Iterate through genotypes to calculate consecutive "1/1" or "1|1"
                for i, genotype in enumerate(genotypes):
                    # Check if the genotype is "1/1" or "1|1"
                    if genotype in counting_genotypes:
                        if previous_row_counts is None: #in case first row
                            current_counts.append(1) #write 1
                        else:
                            # Increment the count from the previous row
                            current_counts.append(previous_row_counts[i] + 1)
                    elif genotype in missing_genotypes:
                        # Missing genotype; carry over the previous count without incrementing
                        if previous_row_counts is None: #in case first row
                            current_counts.append(0) #write 0
                        else:
                            # keep existing count
                            current_counts.append(previous_row_counts[i])
                    else:
                        # Genotype is different; reset the count
                        current_counts.append(0)
                # Append the counts to the new row
                new_row.extend(current_counts)
                result.append(new_row)
                # Save the current row counts as previous row counts for the next iteration
                previous_row_counts = current_counts
    return result

# Function to write a debugging TSV output for individual consecutive data
def write_debug_tsv(consecutive_data, debug_file):
    individual_names = [clean_individual_id(name) for name in consecutive_data[0][2:]]
    positions = [f"{row[0]}:{row[1]}" for row in consecutive_data[1:]]
    
    with open(debug_file, 'w') as tsv:
        # Write the header
        tsv.write("individual\t" + "\t".join(positions) + "\n")
        
        # Write the rows for each individual
        for idx, individual_name in enumerate(individual_names):
            individual_data = [str(row[2 + idx]) for row in consecutive_data[1:]]
            tsv.write(individual_name + "\t" + "\t".join(individual_data) + "\n")

# Combined function to convert consecutive homozygous blocks into a .bed file and generate reports
def convert_and_write_bed_cleaned(consecutive_data, n, n_disrupt, outfile, logfile, debug=False):
    bed_rows = []
    block_frequencies = defaultdict(int)
    individual_block_counts = defaultdict(int)  # To store individual block counts for each strain
    block_length_frequencies = defaultdict(int)  # To store block length distribution (log scale)
    
    # Log-scale bins
    bins = [(1, 10), (10, 100), (100, 1000), (1000, 10000), (10000, 100000), (100000, 1000000), (1000000, 10000000)]
    
    # Process consecutive data by individual
    num_columns = len(consecutive_data[0])  # Number of columns (including chromosome and posistion)
    
    for column_idx in range(2, num_columns):  # Skip chromosome and position columns
        individual_id = clean_individual_id(consecutive_data[0][column_idx])  # Clean individual ID
        individual_block_counts[individual_id] = 0
        current_chromosome = None
        if debug:
          print(f"start individual {individual_id}")
        for row_idx in range(1, len(consecutive_data)):
            chromosome, position = consecutive_data[row_idx][0], int(consecutive_data[row_idx][1])
            consecutive_count = consecutive_data[row_idx][column_idx]
            if debug:
              print(f"start position {position}")
            if current_chromosome is None: #start of new individual
                last_hom_position = None
                last_het_position = None
                block_start = None
                block_end = None
                gap_rows = 0  # Track the number of rows in the disruption
                max_consecutive_count = 0  # Track the number of maximum consecutive count
                current_chromosome = chromosome
                if debug:
                  print(f"start chromosome {current_chromosome}")
                inblock = False
                
            elif current_chromosome != chromosome or row_idx == len(consecutive_data): #when end of a chromosome
              if inblock == True and max_consecutive_count >= n: #end chromosome with a LOH block
                  bed_rows.append([current_chromosome, block_start - 1, block_end, individual_id])
              #reset trackers
              last_hom_position = None
              last_het_position = None
              block_start = None
              block_end = None
              gap_rows = 0
              max_consecutive_count = 0
              inblock = False
              current_chromosome = chromosome
              if debug:
                print(f"end chromosome {current_chromosome}")
                
            if consecutive_count == 0:
              last_het_position = position
              if inblock == True: #count gap
                gap_rows += 1
                if gap_rows == 1: #record potential end block position
                  #block_end = midpoint(last_hom_position, position, round_up=True)
                  block_end = last_hom_position
                if gap_rows >= n_disrupt: # end block
                  inblock = False
                  if max_consecutive_count >= n:
                     bed_rows.append([current_chromosome, block_start - 1, block_end, individual_id])
                  max_consecutive_count = 0
                  
            elif consecutive_count >= 1:
              last_hom_position = position
              inblock = True
              max_consecutive_count = max(max_consecutive_count, consecutive_count)
              block_end = position #in case end with no gap (end of chromosome)
              gap_rows = 0
              if max_consecutive_count == 1: #initial a block, with start position
                if last_hom_position is None:
                  block_start = 1
                else:
                  #block_start = midpoint(last_hom_position, position, round_up=True)
                  block_start = last_hom_position
    for block in bed_rows:
      individual_block_counts[block[3]] += 1
      block_size = int(block[2]) - int(block[1])
      for bin_min, bin_max in bins:
        if bin_min <= block_size < bin_max:
          block_length_frequencies[f"{bin_min}-{bin_max}"] += 1
    for individual_block_number in individual_block_counts.values():
      block_frequencies[individual_block_number] += 1
    
    with open(outfile, 'w') as f:
        for row in bed_rows:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\n")
    
    with open(logfile, 'w') as log:
      # Individual Block Counts
      log.write("Individual Block Counts:\n")
      log.write("Strain\tBlock Count\n")
      print("\nIndividual Block Counts:")
      for individual_id, block_count in individual_block_counts.items():
        log.write(f"{individual_id}\t{block_count}\n")
        print(f"{individual_id}\t{block_count}")
        
      # Block Count Frequencies
      log.write("\nBlock Count Frequencies:\n")
      log.write("Block\tCount\n")
      print("\nBlock Count Frequencies:")
      for count, freq in sorted(block_frequencies.items()):
          log.write(f"{count}\t{freq}\n")
          print(f"{count}\t{freq}")
    
      # Block Length Distribution (Log Scale)
      log.write("\nBlock Length Distribution (Log Scale):\n")
      log.write("Length Range\tFrequency\n")
      print("\nBlock Length Distribution (Log Scale):")
      for length_range, freq in sorted(block_length_frequencies.items()):
          log.write(f"{length_range}\t{freq}\n")
          print(f"{length_range}\t{freq}")


#usage
n_list = [5]  # The threshold for the number of consecutive sites
file_path = 'bwa_haplotypecaller_finalvcf/runs.diploid.vcf.tsv.gz'

# Run the function to process the file with header handling
output_data_with_header = count_consecutive_homozygous_with_header(file_path, ["1/1", "1|1"], [])

for n in n_list:
    outfile = f'LOH_detect/LOH_minSNP-{n}.bed'  # Set the desired output .bed file path
    logfile = f'LOH_detect/LOH_minSNP-{n}.log'  # Set the desired output log file path

    # Convert and write the result to the specified file path and generate the report
    convert_and_write_bed_cleaned(output_data_with_header, n, 30, outfile, logfile)
    
# Run the function to process the file with header handling
output_data_with_header = count_consecutive_homozygous_with_header(file_path, ["0/0", "0|0"], [])

for n in n_list:
    outfile = f'LOH_detect/revLOH_minSNP-{n}.bed'  # Set the desired output .bed file path
    logfile = f'LOH_detect/revLOH_minSNP-{n}.log'  # Set the desired output log file path

    # Convert and write the result to the specified file path and generate the report
    convert_and_write_bed_cleaned(output_data_with_header, n, 30, outfile, logfile)


