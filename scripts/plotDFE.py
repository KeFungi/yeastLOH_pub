import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import re
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=RuntimeWarning)

def parse_s5_params(param_str):
    """
    Parses the detailed parameter string from TableS5_sum_DEF.csv.
    Example: "c=4.60e+08, loc=-1.84e+07, scale=1.84e+07"
    """
    params = {}
    # Use regex to find key=value pairs, handling scientific notation and signs
    matches = re.findall(r'([a-zA-Z_]+)\s*=\s*([-\d.eE+]+)', param_str)
    for key, value in matches:
        params[key] = float(value)
    return params

def parse_s6_params(param_str):
    """
    Parses the parameter string from the TableS6_tail_analysis_summary.csv file.
    Example: "-3.13e-01, -8.26e-03, 1.26e-02"
    """
    # Remove quotes and split by comma, then convert to float
    return [float(p.strip()) for p in param_str.replace('"', '').split(',')]

def plot_whole_dfe(ax, treatment_name, full_df, summary_df, colors):
    """
    Reproduces plots from FigS8. Plots LOH-specific histograms (outline and filled)
    and overlays ALL distribution fits, starring the best one based on AIC.
    **MODIFIED:** Increased font sizes for all text elements.
    """
    # 1. Define LOH conditions and get data subsets
    loh_filters = {
        'LOH=1': (full_df['n_LOH-both'] == 1),
        '1<=LOH<=3': (full_df['n_LOH-both'] >= 1) & (full_df['n_LOH-both'] <= 3)
    }
    data_loh1 = full_df[loh_filters['LOH=1'] & (full_df['treatment'] == treatment_name)]['mean_s'].dropna()
    data_loh3 = full_df[loh_filters['1<=LOH<=3'] & (full_df['treatment'] == treatment_name)]['mean_s'].dropna()

    xmin_data = data_loh3.min()
    xmax_data = data_loh3.max()
    x_range_data = xmax_data - xmin_data
    plot_xmin = xmin_data - 0.05
    plot_xmax = xmax_data + 0.15 * x_range_data
    x_plot = np.linspace(plot_xmin, plot_xmax, 500)

    bins = np.histogram_bin_edges(data_loh3, bins=20)

    # 2. Plot histograms and collect their legend handles
    ax.hist(data_loh3, bins=bins, density=True, label='Strain in MLE (1≤nLOH≤3)',
            histtype='bar', fill=False, edgecolor='gray', lw=1.5)
    ax.hist(data_loh1, bins=bins, density=True, label='Strain nLOH=1',
            color=colors['hist_loh1'], alpha=0.3)
    
    hist_handles, hist_labels = ax.get_legend_handles_labels()
    hist_legend_map = dict(zip(hist_labels, hist_handles))

    final_handles = [hist_legend_map['Strain in MLE (1≤nLOH≤3)'], hist_legend_map['Strain nLOH=1']]
    final_labels = ['Strain in MLE (1≤nLOH≤3)', 'Strain nLOH=1']

    # 3. Find the best model
    treatment_fits = summary_df[summary_df['Treatment'] == treatment_name].copy()
    best_fit_row = treatment_fits.loc[treatment_fits['AIC'].idxmin()]
    best_dist_name = best_fit_row['Distribution']

    dist_map = {
        'Generalized Normal': st.gennorm, 'Weibull': st.weibull_min,
        'Skew-Normal': st.skewnorm, 'Normal': st.norm, 'Gamma': st.gamma,
        'Cauchy': st.cauchy, 'Frechet': st.invweibull,
    }
    
    dist_order = [
        'Generalized Normal', 'Weibull', 'Skew-Normal', 'Normal', 
        'Gamma', 'Cauchy', 'Frechet'
    ]

    # 4. Iterate, plot distributions, and append to legend lists
    pdf_threshold = 0
    for dist_name in dist_order:
        row = treatment_fits[treatment_fits['Distribution'] == dist_name]
        if row.empty:
            continue
        row = row.iloc[0]

        dist = dist_map.get(dist_name)
        if not dist: continue

        params = parse_s5_params(row['Parameters'])
        shape_params = [v for k, v in params.items() if k not in ['loc', 'scale']]
        pdf_vals = dist.pdf(x_plot, *shape_params, loc=params['loc'], scale=params['scale'])
        pdf_vals[pdf_vals < pdf_threshold] = np.nan

        line = ax.plot(x_plot, pdf_vals, lw=2.5, color=colors.get(dist_name, 'black'))[0]
        label = f'* {dist_name}' if dist_name == best_dist_name else f'  {dist_name}'
        
        final_handles.append(line)
        final_labels.append(label)

    # 5. Formatting and creating the legend
    ax.set_title(f'Whole: {treatment_name}', fontsize=22, fontweight='bold')
    ax.set_xlabel('Selection Coefficient (s)', fontsize=18)
    ax.set_ylabel('Density', fontsize=18)
    
    ax.legend(handles=final_handles, labels=final_labels, fontsize=13)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_xlim(plot_xmin, plot_xmax)
    ax.set_ylim(bottom=0)

def plot_tail_dfe(ax, treatment_name, fit_data, tail_summary_data, exp_color, gpd_color, hist_color):
    """
    Reproduces plots from FigS9. Plots LOH-specific tail histograms and overlays
    the Exponential and GPD fits, starring the better model from TableS6.
    **MODIFIED:** Increased font sizes for all text elements.
    """
    # 1. Get tail data subsets
    loh_filters = {
        'LOH=1': (fit_data['n_LOH-both'] == 1),
        '1<=LOH<=3': (fit_data['n_LOH-both'] >= 1) & (fit_data['n_LOH-both'] <= 3)
    }
    data_loh1 = fit_data[loh_filters['LOH=1'] & (fit_data['treatment'] == treatment_name)]['mean_s'].dropna()
    data_loh3 = fit_data[loh_filters['1<=LOH<=3'] & (fit_data['treatment'] == treatment_name)]['mean_s'].dropna()
    tail_loh1 = np.sort(data_loh1)[-20:]
    tail_loh3 = np.sort(data_loh3)[-20:]

    xmin_data = tail_loh3.min()
    xmax_data = tail_loh3.max()
    x_range_data = xmax_data - xmin_data
    plot_xmin = xmin_data
    plot_xmax = xmax_data + 0.20 * x_range_data
    x_plot = np.linspace(plot_xmin, plot_xmax, 500)

    bins = np.histogram_bin_edges(tail_loh3, bins=6)

    # 2. Plot histograms
    ax.hist(tail_loh3, bins=bins, density=True, label='Strain in MLE (1≤nLOH≤3)',
            histtype='bar', fill=False, edgecolor='gray', lw=1.5)
    ax.hist(tail_loh1, bins=bins, density=True, label='Strain nLOH=1',
            color=hist_color, alpha=0.3)

    # 3. Get parameters
    row = tail_summary_data[tail_summary_data['Treatment'] == treatment_name].iloc[0]
    better_model = row['Better_Model']
    exp_params_list = parse_s6_params(row['Exp_Params (loc, scale)'])
    gpd_params_list = parse_s6_params(row['GPD_Params (c, loc, scale)'])
    exp_params = {'loc': exp_params_list[0], 'scale': exp_params_list[1]}
    gpd_params = {'c': gpd_params_list[0], 'loc': gpd_params_list[1], 'scale': gpd_params_list[2]}

    # 4. Plot PDFs
    pdf_threshold = 1e-10
    exp_pdf = st.expon.pdf(x_plot, **exp_params)
    gpd_pdf = st.genpareto.pdf(x_plot, **gpd_params)
    exp_pdf[exp_pdf < pdf_threshold] = np.nan
    gpd_pdf[gpd_pdf < pdf_threshold] = np.nan

    exp_label_base = 'Exponential'
    gpd_label_base = 'GPD'
    exp_label = f'* {exp_label_base}' if 'Exponential' in better_model else f'  {exp_label_base}'
    gpd_label = f'* {gpd_label_base}' if 'GPD' in better_model else f'  {gpd_label_base}'

    ax.plot(x_plot, exp_pdf, color=exp_color, lw=2.5, label=exp_label)
    ax.plot(x_plot, gpd_pdf, color=gpd_color, lw=2.5, label=gpd_label)

    # 5. Formatting
    ax.set_title(f'Tail: {treatment_name}', fontsize=22, fontweight='bold')
    ax.set_xlabel('Selection Coefficient (s)', fontsize=18)
    ax.set_ylabel('Density', fontsize=18)
    
    handles, labels = ax.get_legend_handles_labels()
    
    order = ['Strain in MLE (1≤nLOH≤3)', 'Strain nLOH=1']
    if 'Exponential' in better_model:
        order.extend(['* Exponential', '  GPD'])
    else:
        order.extend(['  Exponential', '* GPD'])

    legend_map = dict(zip(labels, handles))
    sorted_handles = [legend_map[lbl] for lbl in order if lbl in legend_map]
    sorted_labels = [lbl for lbl in order if lbl in legend_map]
    
    ax.legend(handles=sorted_handles, labels=sorted_labels, fontsize=14, loc='upper right')
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.set_xlim(plot_xmin, plot_xmax)
    ax.set_ylim(bottom=0)

def main():
    """ Main function to load data, generate all plots, and save to a single PDF. """
    # --- 1. Setup ---
    khroma_muted = {
        'hist_loh1_unified': '#E69F00',
        'Weibull': '#332288', 'Skew-Normal': '#AA4499', 'Generalized Normal': '#44AA99',
        'Normal': '#117733', 'Frechet': '#999933', 'Gamma': '#DDCC77', 'Cauchy': '#882255',
        'exp': '#CC6677', 'gpd': '#332288'
    }

    try:
        fit_tb = pd.read_csv('tables/fit_tb.csv')
        sum_def = pd.read_csv('TableS5_sum_DEF.csv')
        tail_summary = pd.read_csv('TableS6_tail_analysis_summary.csv')
    except FileNotFoundError as e:
        print(f"Error: Could not find input file - {e.filename}")
        return

    output_pdf_filename = 'plots/reproduced_dfe_plots_final.pdf'
    
    # --- 2. Create Plots ---
    with PdfPages(output_pdf_filename) as pdf:
        print(f"Generating 9 plots and saving to '{output_pdf_filename}'...")
        fig, axes = plt.subplots(3, 3, figsize=(20, 20))
        axes = axes.flatten()

        plot_colors = {k: v for k, v in khroma_muted.items()}
        plot_colors['hist_loh1'] = khroma_muted['hist_loh1_unified']

        plot_whole_dfe(axes[0], 'Caffeine', fit_tb, sum_def, plot_colors)
        plot_whole_dfe(axes[1], 'SC', fit_tb, sum_def, plot_colors)

        tail_treatments = tail_summary['Treatment'].unique()
        for i, treatment in enumerate(tail_treatments):
            if i + 2 < len(axes):
                plot_tail_dfe(axes[i+2], treatment, fit_tb, tail_summary, 
                              exp_color=khroma_muted['exp'], 
                              gpd_color=khroma_muted['gpd'],
                              hist_color=khroma_muted['hist_loh1_unified'])
        
        for i in range(len(tail_treatments) + 2, len(axes)):
            axes[i].set_visible(False)

        fig.tight_layout(rect=[0, 0.03, 1, 0.96])
        pdf.savefig(fig)
        plt.close(fig)
        print("PDF generation complete.")

if __name__ == '__main__':
    main()

