import pandas as pd
import numpy as np
import scipy.stats as st
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import warnings
import concurrent.futures
import os
from matplotlib.backends.backend_pdf import PdfPages

# --- Helper Function for Robust Initial Guesses ---

def get_robust_initial_guess(name, dist, data):
    """
    Generates a more stable initial guess for notoriously difficult distributions.
    """
    try:
        if name == 'Cauchy':
            # For Cauchy, median and IQR are robust estimators for loc and scale
            loc = np.median(data)
            scale = (np.quantile(data, 0.75) - np.quantile(data, 0.25)) / 2
            return (loc, scale)
        elif name == 'Gamma':
            # For Gamma, the data must be > loc. We estimate loc first.
            loc = np.min(data) - 1e-4 # Estimate loc as just below the minimum
            shifted_data = data - loc
            # Fit shape (a) and scale on the shifted data, with loc fixed at 0
            a, _, scale = st.gamma.fit(shifted_data, floc=0)
            return (a, loc, scale)
        else:
            # For other distributions, the standard fit is usually sufficient
            return dist.fit(data)
    except Exception as e:
        print(f"  -> Initial guess failed for {name}: {e}. Using default guess.")
        # Fallback to a generic guess if everything fails
        return [1.0] * dist.numargs + [np.mean(data), np.std(data)]

# --- Monte Carlo Likelihood Function for Sum Distributions ---

def monte_carlo_sum_likelihood(params, n, x_vals, dist, num_samples=1000000, bin_size=100):
    """
    Numerically estimates the PDF for the sum of 'n' random variables from a given distribution.
    """
    shape_params = params[:-2]
    loc = params[-2]
    scale = params[-1]

    samples = dist.rvs(*shape_params, loc=loc, scale=scale, size=(n, int(num_samples)))
    sum_samples = np.sum(samples, axis=0)
    sum_samples = sum_samples[np.isfinite(sum_samples)]
    if len(sum_samples) == 0:
        return np.full_like(x_vals, 1e-10, dtype=float)

    combined_data = np.concatenate((sum_samples, x_vals))
    min_val = np.min(combined_data)
    max_val = np.max(combined_data)
    
    hist, bin_edges = np.histogram(sum_samples, range=(min_val, max_val), bins=bin_size, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    likelihood = np.interp(x_vals, bin_centers, hist)
    likelihood[likelihood == 0] = 1e-10
    return likelihood

# --- MLE Function Using the Monte Carlo Likelihood ---

def fit_sum_distribution(data, event_counts, dist, initial_params):
    """
    Finds the MLE parameters for a distribution where data points are sums of 'n' events.
    """
    def negative_log_likelihood(params, data, event_counts, dist):
        log_likelihood = 0
        
        precomputed_pdfs = {}
        for n in np.unique(event_counts):
            if n == 0: continue
            
            x_vals = np.unique(data[event_counts == n])
            l_vals = monte_carlo_sum_likelihood(params, int(n), x_vals, dist)
            precomputed_pdfs[n] = {x: l for x, l in zip(x_vals, l_vals)}

        for d, n in zip(data, event_counts):
            if n == 0: continue
            likelihood = precomputed_pdfs[n][d]
            log_likelihood += np.log(likelihood)
            
        if not np.isfinite(log_likelihood):
            return np.inf

        return -log_likelihood
    
    # Define parameter bounds to aid optimization
    bounds = []
    # Shape parameters (often > 0, but skew-normal can be negative)
    if dist.name == 'skewnorm':
        bounds.append((None, None))
    else:
        for _ in range(dist.numargs):
            bounds.append((1e-5, None)) 
    # Location parameter (can be anything)
    bounds.append((None, None))
    # Scale parameter (must be > 0)
    bounds.append((1e-5, None))

    result = minimize(
        negative_log_likelihood,
        x0=initial_params,
        args=(data, event_counts, dist),
        method='L-BFGS-B',
        bounds=bounds
    )

    if result.success and np.isfinite(result.fun):
        final_nll = result.fun
        k = len(result.x)
        aic = 2 * k + 2 * final_nll
        return {'params': result.x, 'aic': aic}
    else:
        return {'params': None, 'aic': np.inf}

# --- MODIFIED: Worker Function to try until 3 successes or 10 attempts ---

def worker_fit_distribution(args):
    """
    A wrapper function for a single fitting task. It tries until it gets
    3 successful fits or reaches 10 total attempts, then returns the best result.
    """
    name, dist, data, event_counts = args
    print(f"Fitting {name} distribution on process {os.getpid()}...")
    
    MAX_ATTEMPTS = 10
    REQUIRED_SUCCESSES = 3
    
    successful_fits = []
    attempt_count = 0

    while len(successful_fits) < REQUIRED_SUCCESSES and attempt_count < MAX_ATTEMPTS:
        attempt_count += 1
        
        # Generate a new initial guess for each attempt
        initial_params = get_robust_initial_guess(name, dist, data)
        
        current_result = fit_sum_distribution(data, event_counts, dist, initial_params)

        # --- Print intermediate results for each attempt ---
        print(f"  -> Attempt {attempt_count}/{MAX_ATTEMPTS} for {name}:")
        if current_result['params'] is not None:
            param_str = ", ".join([f"{p:.4f}" for p in current_result['params']])
            print(f"     SUCCESS | AIC: {current_result['aic']:.2f}, Params: [{param_str}]")
            successful_fits.append(current_result)
        else:
            print("     Fit FAILED.")
    
    # --- Select the best fit from all successful attempts ---
    if not successful_fits:
        print(f"  -> All {attempt_count} fitting attempts for {name} FAILED.")
        best_result = {'params': None, 'aic': np.inf}
    else:
        # Find the result with the minimum 'aic' value (highest likelihood)
        best_result = min(successful_fits, key=lambda x: x['aic'])
        print(f"  -> Best fit for {name} selected from {len(successful_fits)} successful attempts (Best AIC: {best_result['aic']:.2f}).")
    
    return name, best_result

# --- Main Analysis Script ---

def main():
    warnings.filterwarnings("ignore")
    
    # Create directories for output if they don't exist
    os.makedirs('tables', exist_ok=True)
    os.makedirs('plots', exist_ok=True)

    try:
        df = pd.read_csv('tables/fit_tb.csv')
    except FileNotFoundError:
        print("Error: 'tables/fit_tb.csv' not found.")
        return

    df['datapoint'] = df['mean_s']
    df['event_count'] = df['n_LOH-both']
    
    distributions_to_fit = {
        'Normal': st.norm,
        'Weibull': st.weibull_min,
        'Frechet': st.invweibull,
        'Gamma': st.gamma,
        'Cauchy': st.cauchy,
        'Generalized Normal': st.gennorm,
        'Skew-Normal': st.skewnorm
    }

    treatments_to_analyze = ['Caffeine', 'SC']
    df_filtered = df[df['treatment'].isin(treatments_to_analyze) & (df['event_count'] >= 1) & (df['event_count'] <= 3)].copy()

    all_results = {}
    for treatment in treatments_to_analyze:
        print(f"\n--- Processing Treatment: {treatment} ---")
        subset_df = df_filtered[df_filtered['treatment'] == treatment]
        data = subset_df['datapoint'].dropna().to_numpy()
        event_counts = subset_df['event_count'].dropna().to_numpy()

        if len(data) == 0:
            print("No valid data for this treatment after filtering.")
            continue
        
        tasks = [(name, dist, data, event_counts) for name, dist in distributions_to_fit.items()]
        
        treatment_results = {}
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results_iterator = executor.map(worker_fit_distribution, tasks)
            for name, result in results_iterator:
                treatment_results[name] = result
        
        all_results[treatment] = treatment_results

    # Print summary
    summary_data = []
    for treatment, treatment_results in all_results.items():
        sorted_results = sorted(treatment_results.items(), key=lambda item: item[1]['aic'])
        for name, res in sorted_results:
            if res['params'] is not None:
                dist = distributions_to_fit[name]
                param_names = []
                if dist.shapes:
                    param_names.extend(dist.shapes.split(','))
                param_names.extend(['loc', 'scale'])
                
                param_str_parts = []
                for p_name, p_val in zip(param_names, res['params']):
                    param_str_parts.append(f"{p_name}={p_val}")
                param_str = ', '.join(param_str_parts)
            else:
                param_str = 'Fit Failed'

            summary_data.append({
                'Treatment': treatment, 'Distribution': name,
                'AIC': f"{res['aic']:.2f}" if res['aic'] != np.inf else 'N/A',
                'Parameters': param_str
            })
    summary_df = pd.DataFrame(summary_data)
    print("\n--- Summary of Best Sum Distribution Fits (Lower AIC is better) ---")
    print(summary_df.to_string())
    print("--------------------------------------------------------------------")

    summary_df.to_csv('tables/sum_DEF.csv', index=False)
    print("\nResults table saved to 'tables/sum_DEF.csv'")


    # --- Plotting Section ---
    plt.rcParams.update({'font.size': 14})
    
    pdf_filename = 'plots/sum_DEF_combined.pdf'
    with PdfPages(pdf_filename) as pdf:
        # Create a single figure with two subplots vertically stacked
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 16))
        
        for i, treatment in enumerate(treatments_to_analyze):
            ax = axes[i] # Get the current subplot (ax)
            
            treatment_results = all_results.get(treatment)
            if not treatment_results:
                ax.text(0.5, 0.5, f'No data for {treatment}', ha='center', va='center')
                continue

            all_events_df = df[(df['treatment'] == treatment) & (df['event_count'] >= 1) & (df['event_count'] <= 3)]
            all_events_data = all_events_df['datapoint'].dropna().to_numpy()

            subset_df = df_filtered[df_filtered['treatment'] == treatment]
            data_n1 = subset_df[subset_df['event_count'] == 1]['datapoint'].dropna().to_numpy()

            if len(all_events_data) == 0:
                print(f"Skipping plot for {treatment}: No data to display.")
                ax.text(0.5, 0.5, f'No data for {treatment}', ha='center', va='center')
                continue
            
            plot_max = np.quantile(all_events_data, 0.999) * 1.1
            plot_min = np.quantile(all_events_data, 0.001) * 1.1
            
            ax.hist(all_events_data, bins=25, density=True, label='Data Points in Analysis (LOH ≤ 3)', 
                     range=(plot_min, plot_max), histtype='bar', facecolor='none', 
                     edgecolor='black', linewidth=1.2, zorder=2)
            
            if len(data_n1) > 0:
                ax.hist(data_n1, bins=25, density=True, alpha=0.6, 
                         label='Data Points LOH = 1', range=(plot_min, plot_max), zorder=1)

            sorted_results = sorted(treatment_results.items(), key=lambda item: item[1]['aic'])
            
            for idx, (name, res) in enumerate(sorted_results):
                if res['params'] is not None:
                    dist = distributions_to_fit[name]
                    label = name
                    
                    if idx == 0:
                        label = f"* {label}"

                    x = np.linspace(plot_min, plot_max, 1000)

                    # --- Correctly unpack shape, loc, and scale parameters ---
                    # The number of shape parameters is given by dist.numargs
                    num_shape_params = dist.numargs
                    shape_params = res['params'][:num_shape_params]
                    loc_param = res['params'][num_shape_params]
                    scale_param = res['params'][num_shape_params + 1]
                    
                    pdf_vals = dist.pdf(x, *shape_params, loc=loc_param, scale=scale_param)
                    
                    ax.plot(x, pdf_vals, label=label, linewidth=2, zorder=3)

            ax.set_title(f'{treatment}',  fontsize=25, fontweight='bold')
            ax.set_xlabel('s', fontsize=20)
            ax.set_ylabel('Density', fontsize=20)
            ax.tick_params(axis='both', which='major', labelsize=16)
            ax.legend(fontsize=12, loc='upper left')
            ax.grid(True)
            ax.set_xlim(plot_min, plot_max)
            ax.set_ylim(bottom=0)

        fig.tight_layout(pad=4.0)
        pdf.savefig(fig)
        print(f"\nCombined plot saved to '{pdf_filename}'")
        
        for i, ax in enumerate(axes.flatten()):
            treatment = treatments_to_analyze[i]
            extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
            png_filename = f'plots/sum_DFE_{treatment}.png'
            fig.savefig(png_filename, bbox_inches=extent)
            print(f"Generated individual plot: {png_filename}")

        plt.close(fig)

if __name__ == "__main__":
    main()

