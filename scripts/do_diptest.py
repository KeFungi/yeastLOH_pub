import pandas as pd
import diptest
import warnings

def perform_diptest_on_strains():
    """
    Loads data, filters out control strains, and performs a dip test for bimodality
    on specified treatments, then prints the results in a table.
    """
    warnings.filterwarnings("ignore", category=FutureWarning)

    # 1. Load fit_tb.csv
    try:
        fit_tb = pd.read_csv('tables/fit_tb.csv')
    except FileNotFoundError:
        print("Error: 'tables/fit_tb.csv' not found. Please ensure the file is in the current directory.")
        return

    # 2. Use the same unique_treatment list
    unique_treatments = ["Acidic", "Caffeine", "H2O2", "SC", "worm-20C", "YPD-30C", "YPD-37C", "(average)"]
    
    # 3. Filter data for strains that DO NOT have the 'CNTRL-' prefix
    # The `na=False` ensures that any NaN values in the 'strain' column are treated as not matching.
    filtered_tb = fit_tb[~fit_tb['strain'].str.startswith('CNTRL-', na=False)].copy()

    results_list = []

    print("--- Performing Dip Test on Strains (excluding 'CNTRL-') ---")
    
    # 4. Do dip test for each treatment
    for treatment in unique_treatments:
        
        treatment_data = filtered_tb[filtered_tb['treatment'] == treatment]['mean_s'].dropna()
        
        sample_size = len(treatment_data)
        
        if sample_size > 3:  # Diptest requires at least 4 data points
            dip_statistic, p_value = diptest.diptest(treatment_data.to_numpy())
            
            # 5. Record results
            results_list.append({
                'Treatment': treatment,
                'Sample Size': sample_size,
                'Dip Statistic': f"{dip_statistic:.4f}",
                'P-Value': f"{p_value:.3e}"
            })
        else:
            results_list.append({
                'Treatment': treatment,
                'Sample Size': sample_size,
                'Dip Statistic': 'N/A',
                'P-Value': 'N/A'
            })

    # Convert the list of results into a pandas DataFrame for nice printing
    results_df = pd.DataFrame(results_list)
    
    print("\nDiptest Results:")
    print(results_df.to_string(index=False))
    
    results_df.to_csv('tables/diptest.csv', index=False)

if __name__ == "__main__":
    perform_diptest_on_strains()
