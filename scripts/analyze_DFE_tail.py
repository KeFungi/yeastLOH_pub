import pandas as pd
import numpy as np
import scipy.stats as st
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import warnings
import concurrent.futures
import os
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import find_peaks

# --- Analysis Configuration ---
NUM_TAIL_POINTS = 20
MAX_CPU_CORES = 4

# --- Helper Functions ---

def likelihood_ratio_test(ll_null, ll_alternative, df):
    """Performs a Likelihood-Ratio Test."""
    lrt_stat = 2 * (ll_alternative - ll_null)
    p_value = st.chi2.sf(lrt_stat, df)
    return lrt_stat, p_value

def get_gpd_type(c):
    """Classifies the GPD type based on the shape parameter 'c'."""
    if c < 0:
        return "GPD (Weibull)"
    elif c > 0:
        return "GPD (Frechet)"
    else: # c == 0
        return "GPD (Exponential)"

# --- Trapezoidal Integration for OPTIMIZATION & LRT CALCULATION ---

class TrapezoidalConvolution:
    """Manages recursive convolution calculation via fast `trapz` integration."""
    def __init__(self, dist, params, data_range, linspace_params):
        self.dist, self._pdf_cache = dist, {}
        if self.dist.name == 'genpareto':
            self.shape_params, self.loc, self.scale = params[:1], params[1], params[2]
        else:
            self.shape_params, self.loc, self.scale = [], params[0], params[1]
        
        # Use provided linspace parameters for the integration grid
        self.x_grid = np.linspace(*linspace_params)


    def get_pdf(self, n):
        if n in self._pdf_cache: return self._pdf_cache[n]
        if n == 1:
            pdf_values = self.dist.pdf(self.x_grid, *self.shape_params, loc=self.loc, scale=self.scale)
            self._pdf_cache[1] = lambda x: np.interp(x, self.x_grid, pdf_values, left=0, right=0)
            return self._pdf_cache[1]
            
        pdf_n_minus_1_func = self.get_pdf(n - 1)
        pdf_1_func = self.get_pdf(1)
        pdf_n_minus_1_vals = pdf_n_minus_1_func(self.x_grid)
        pdf_1_vals = pdf_1_func(self.x_grid)

        def convoluted_pdf(z):
            pdf_1_shifted_vals = np.interp(z - self.x_grid, self.x_grid, pdf_1_vals, left=0, right=0)
            integrand = pdf_n_minus_1_vals * pdf_1_shifted_vals
            return np.trapz(integrand, self.x_grid)
            
        self._pdf_cache[n] = convoluted_pdf
        return convoluted_pdf

# --- Log-likelihood calculation (penalty removed) ---
def calculate_loglik_trapz(params, data, event_counts, dist, linspace_params):
    """Calculates the final log-likelihood using trapezoidal integration."""
    log_likelihood = 0
    conv_manager = TrapezoidalConvolution(dist, params, data, linspace_params)
    pdf_functions = {n: conv_manager.get_pdf(n) for n in np.unique(event_counts) if n > 0}
    for i, (d, n) in enumerate(zip(data, event_counts)):
        if n == 0: continue
        pdf_func = pdf_functions.get(n)
        if pdf_func:
            likelihood = pdf_func(d)
            if likelihood > 1e-100: 
                log_likelihood += np.log(likelihood)
            else: 
                log_likelihood -= 100
    return log_likelihood if np.isfinite(log_likelihood) else -np.inf

# --- MODIFIED: Optimization function now accepts a hard location boundary ---
def fit_sum_distribution_trapz(tail_data, tail_event_counts, dist, initial_params=None, max_loc=None, linspace_params=None):
    """Finds MLE parameters and log-likelihood using TRAPEZOIDAL integration for the optimization."""
    def negative_log_likelihood_trapz(params, data, event_counts, dist, linspace_params):
        loglik = calculate_loglik_trapz(params, data, event_counts, dist, linspace_params)
        print(f"    Trying params: {['{:.4e}'.format(p) for p in params]}, LogLik: {loglik:.4f}", flush=True)
        return -loglik if np.isfinite(loglik) else np.inf

    if dist.name == 'genpareto':
        if initial_params is None:
            initial_loc = np.mean(tail_data)
            # --- FIX: Clip the initial guess to be within the bounds ---
            if max_loc is not None:
                initial_loc = min(initial_loc, max_loc)
            initial_params = [0, initial_loc, np.std(tail_data)]
        bounds = [(-1, None), (None, max_loc), (1e-5, None)]
    else: # Exponential
        if initial_params is None:
            initial_loc = np.mean(tail_data)
            # --- FIX: Clip the initial guess to be within the bounds ---
            if max_loc is not None:
                initial_loc = min(initial_loc, max_loc)
            initial_params = [initial_loc, np.std(tail_data)]
        bounds = [(None, max_loc), (1e-5, None)]

    # --- Final check to ensure provided initial_params also respect bounds ---
    if max_loc is not None:
        loc_index = -2 if dist.name != 'genpareto' else 1
        if initial_params[loc_index] > max_loc:
            initial_params[loc_index] = max_loc

    result = minimize(
        negative_log_likelihood_trapz, x0=initial_params,
        args=(tail_data, tail_event_counts, dist, linspace_params),
        method='L-BFGS-B',
        bounds=bounds
    )

    if result.success and np.isfinite(result.fun):
        return {'params': result.x, 'loglik': -result.fun}
    else:
        return {'params': None, 'loglik': -np.inf}

# --- MODIFIED Worker Function ---
def worker_process_treatment(args):
    """
    Processes a single treatment: fits Exponential and GPD models independently,
    enforcing a hard bound on the location parameter. Tries different integration
    settings if the first one fails.
    """
    treatment, tail_data, tail_event_counts = args
    
    # Define different integration grid settings to try
    linspace_opts = [
        (-1, 1, 1000000),       # First attempt
        (-0.5, 0.5, 10000000)   # Second attempt if first fails
    ]

    # Initialize results
    exp_fit_result = {'params': None, 'loglik': -np.inf}
    gpd_fit_result = {'params': None, 'loglik': -np.inf}
    final_results = {}

    # --- Determine the maximum possible location from n=1 data for the hard bound ---
    max_loc_bound = None
    if 1 in tail_event_counts:
        n1_data = tail_data[tail_event_counts == 1]
        if len(n1_data) > 0:
            max_loc_bound = np.min(n1_data)
            print(f"  -> For {treatment}, location parameter will be bounded at max: {max_loc_bound:.4e}", flush=True)

    # Loop through each set of linspace parameters as an "attempt"
    for i, ls_params in enumerate(linspace_opts):
        print(f"\n--- Starting Attempt {i+1}/{len(linspace_opts)} for {treatment} ---", flush=True)

        # --- FIT EXPONENTIAL ---
        print(f"  -> Optimizing Exponential...", flush=True)
        current_exp_fit = fit_sum_distribution_trapz(
            tail_data, tail_event_counts, st.expon,
            max_loc=max_loc_bound, linspace_params=ls_params
        )

        if current_exp_fit['params'] is None:
            print(f"  -> Exponential fit FAILED on attempt {i+1}. Moving to next attempt.", flush=True)
            continue

        # --- FIT GPD ---
        print(f"  -> Optimizing GPD...", flush=True)
        current_gpd_fit = fit_sum_distribution_trapz(
            tail_data, tail_event_counts, st.genpareto,
            max_loc=max_loc_bound, linspace_params=ls_params
        )

        if current_gpd_fit['params'] is None:
            print(f"  -> GPD fit FAILED on attempt {i+1}. Moving to next attempt.", flush=True)
            continue

        # --- SUCCESS ---
        print(f"  -> Both models SUCCEEDED on attempt {i+1}.", flush=True)
        exp_fit_result = current_exp_fit
        gpd_fit_result = current_gpd_fit
        break # Exit the loop since we have a successful pair of fits.

    # Package the final results
    final_results['Exponential'] = exp_fit_result
    final_results['Generalized Pareto'] = gpd_fit_result

    if exp_fit_result['params'] is None or gpd_fit_result['params'] is None:
         print(f"  -> Fitting FAILED for {treatment} after all attempts.", flush=True)
    
    return treatment, final_results

# --- Main Analysis Script ---

def main():
    warnings.filterwarnings("ignore", category=UserWarning) 
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    try:
        df = pd.read_csv('tables/fit_tb.csv')
    except FileNotFoundError:
        print("Error: 'tables/fit_tb.csv' not found.", flush=True)
        return

    df['datapoint'] = df['mean_s']
    df['event_count'] = df['n_LOH-both']

    unique_treatments = ["Acidic", "Caffeine", "H2O2", "SC", "worm-20C", "YPD-30C", "YPD-37C"]
    print(f"--- Preparing Tasks for Tail Analysis ---", flush=True)
    
    tasks, skipped_treatments = [], []
    for treatment in unique_treatments:
        # Filter for 1 <= n <= 3 BEFORE selecting the tail
        filtered_df = df[
            (df['treatment'] == treatment) & 
            (df['event_count'] >= 1) & 
            (df['event_count'] <= 3)
        ].dropna(subset=['datapoint', 'event_count'])
        
        # Check if there are enough points in this pre-filtered set
        if len(filtered_df) < NUM_TAIL_POINTS:
            print(f"  -> Skipping {treatment}: Insufficient data points ({len(filtered_df)}) with 1<=n<=3 to select a tail of {NUM_TAIL_POINTS}.", flush=True)
            skipped_treatments.append(treatment)
            continue
        
        # Now select the tail from the pre-filtered data
        analysis_df = filtered_df.nlargest(NUM_TAIL_POINTS, 'datapoint')

        if analysis_df.empty:
            print(f"  -> Skipping {treatment}: No data points available after filtering and tail selection.", flush=True)
            skipped_treatments.append(treatment)
            continue

        tail_data = analysis_df['datapoint'].to_numpy()
        tail_event_counts = analysis_df['event_count'].to_numpy().astype(int)
        
        tasks.append((treatment, tail_data, tail_event_counts))

    all_results = {t: {} for t in unique_treatments if t not in skipped_treatments}
    
    print("\n--- Starting Parallel Fitting Process (Trapezoidal Optimization) ---", flush=True)
    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_CPU_CORES) as executor:
            for treatment, treatment_results in executor.map(worker_process_treatment, tasks):
                all_results[treatment] = treatment_results
    except Exception as e:
        print(f"\n--- ERROR during parallel processing: {e} ---", flush=True)

    print("\n--- Intermediate Results ---", flush=True)
    for treatment in unique_treatments:
        if treatment in skipped_treatments: continue
        results = all_results.get(treatment, {})
        if 'Exponential' in results and 'Generalized Pareto' in results and \
           results['Exponential']['params'] is not None and results['Generalized Pareto']['params'] is not None:
            exp_res, gp_res = results['Exponential'], results['Generalized Pareto']
            lrt_stat, p_val = likelihood_ratio_test(exp_res['loglik'], gp_res['loglik'], df=1)
            print(f"  -> {treatment}: LRT Stat={lrt_stat:.2f}, P-Value={p_val:.3f}", flush=True)
        else:
            print(f"  -> Fit failed for one or more models for {treatment}.", flush=True)

    summary_data = []
    for treatment, results in all_results.items():
        fit_failed = not results or any(res.get('params') is None for res in results.values())
        
        filtered_df = df[
            (df['treatment'] == treatment) & 
            (df['event_count'] >= 1) & 
            (df['event_count'] <= 3)
        ].dropna(subset=['datapoint', 'event_count'])
        analysis_df = filtered_df.nlargest(NUM_TAIL_POINTS, 'datapoint')

        if fit_failed:
            summary_data.append({'Treatment': treatment, 'Better_Model': 'Fit Failed'})
            continue

        exp_result = results['Exponential']
        gp_result = results['Generalized Pareto']
        lrt_stat, p_value = likelihood_ratio_test(exp_result['loglik'], gp_result['loglik'], df=1)
        gpd_shape_c = gp_result['params'][0]
        gpd_type = get_gpd_type(gpd_shape_c)
        better_model = gpd_type if p_value < 0.05 else "Exponential"
        
        summary_data.append({
            'Treatment': treatment, 'Threshold': analysis_df['datapoint'].min(),
            'Tail_Points': len(analysis_df),
            'Exp_Params (loc, scale)': f"{exp_result['params'][0]}, {exp_result['params'][1]}",
            'GPD_Type': gpd_type,
            'GPD_Params (c, loc, scale)': f"{gp_result['params'][0]}, {gp_result['params'][1]}, {gp_result['params'][2]}",
            'Exp_LogLik': exp_result['loglik'], 'GPD_LogLik': gp_result['loglik'],
            'P_Value': p_value, 'Better_Model': better_model
        })
        
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        print("\n--- Summary of Tail Distribution Fits ---", flush=True)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000):
            print(summary_df, flush=True)
        summary_df.to_csv("tables/tail_analysis_summary.csv", index=False)
        print("\nTail analysis summary saved to 'tail_analysis_summary.csv'", flush=True)

    pdf_filename = 'plots/tail_analysis_plots.pdf'
    with PdfPages(pdf_filename) as pdf:
        num_treatments_to_plot = len(all_results)
        
        ncols = 2
        nrows = (num_treatments_to_plot + ncols - 1) // ncols
        
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 7 * nrows), constrained_layout=True)
        axes = axes.flatten()
        
        for i, (treatment, results) in enumerate(all_results.items()):
            ax = axes[i]
            
            filtered_df = df[
                (df['treatment'] == treatment) & 
                (df['event_count'] >= 1) & 
                (df['event_count'] <= 3)
            ].dropna(subset=['datapoint', 'event_count'])
            analysis_df = filtered_df.nlargest(NUM_TAIL_POINTS, 'datapoint')
            tail_data_for_analysis = analysis_df['datapoint'].to_numpy()
            
            tail_data_n1 = analysis_df[analysis_df['event_count'] == 1]['datapoint'].to_numpy()
            
            data_min = tail_data_for_analysis.min()
            data_max = tail_data_for_analysis.max()
            common_bins = np.linspace(data_min, data_max, 7)

            ax.hist(tail_data_for_analysis, bins=common_bins, density=True, label='Tail Points in Analysis (LOH ≤ 3)', histtype='bar', facecolor='none', edgecolor='black', linewidth=1.2, zorder=2)
            if len(tail_data_n1) > 0: ax.hist(tail_data_n1, bins=common_bins, density=True, alpha=0.6, label='Tail Points LOH = 1', zorder=1)
            
            ax.set_xlabel("s", fontsize=20)
            ax.set_ylabel("Density", fontsize=20)
            ax.tick_params(axis='both', which='major', labelsize=16)
            ax.grid(True)
            
            fit_succeeded = results and all(res.get('params') is not None for res in results.values())
            
            if fit_succeeded:
                summary_row = summary_df[summary_df['Treatment'] == treatment]
                better_model_name = ""
                if not summary_row.empty:
                    better_model_name = summary_row['Better_Model'].iloc[0]

                exp_label = 'Exponential PDF'
                gpd_label = 'GPD PDF'
                if 'Exponential' in better_model_name:
                    exp_label = '* ' + exp_label
                elif 'GPD' in better_model_name:
                    gpd_label = '* ' + gpd_label
                
                title = f"{treatment}"
                
                x_plot = np.linspace(data_min, data_max, 1000)
                
                # --- MODIFICATION START ---
                # Calculate PDFs
                exp_pdf = st.expon.pdf(x_plot, *results['Exponential']['params'])
                gp_pdf = st.genpareto.pdf(x_plot, c=results['Generalized Pareto']['params'][0], loc=results['Generalized Pareto']['params'][1], scale=results['Generalized Pareto']['params'][2])
                
                # Find the first index where the PDF is positive to avoid plotting the zero-line
                start_index_exp = np.argmax(exp_pdf > 0)
                start_index_gp = np.argmax(gp_pdf > 0)

                # Plot the PDFs starting from their first non-zero value
                ax.plot(x_plot[start_index_exp:], exp_pdf[start_index_exp:], 'r-', lw=2, label=exp_label)
                ax.plot(x_plot[start_index_gp:], gp_pdf[start_index_gp:], 'g-', lw=2, label=gpd_label)
                # --- MODIFICATION END ---
                
                ax.set_title(title, fontsize=25, fontweight='bold')
                ax.legend(loc='upper right', fontsize=16)
            else:
                ax.set_title(f"{treatment}")
                ax.legend(loc='upper right', fontsize=16)

        for i in range(num_treatments_to_plot, len(axes)):
            axes[i].set_visible(False)
        
        pdf.savefig(fig)
        print(f"\nCombined plot saved to '{pdf_filename}'", flush=True)
        plt.close(fig)

if __name__ == "__main__":
    main()

