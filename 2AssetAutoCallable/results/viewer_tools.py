import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import glob
import os
import ipywidgets as widgets
from IPython.display import display, clear_output

# Function to extract parameters from filename
def extract_params(filename):
    match = re.search(r'results_(.*)_run_([\d\.]+)_(\d+)\.csv', filename)
    if match:
        ladder = match.group(1)
        smoothing = float(match.group(2))
        paths = int(match.group(3))
        return smoothing, paths, ladder
    else:
        return None, None, None

# Load all CSV files and organize by smoothing and paths
def load_all_data(file_pattern='sensitivity_results_*.csv'):
    data_by_params = {}
    file_list = glob.glob(file_pattern)
    
    for file_path in file_list:
        filename = os.path.basename(file_path)
        smoothing, paths, ladder = extract_params(filename)
        
        if smoothing is not None and paths is not None and ladder is not None:
            df = pd.read_csv(file_path)
            data_by_params[(smoothing, paths, ladder)] = df
            
    return data_by_params

# Create widgets for interactive exploration
def create_interactive_explorer(ladders, smoothing_params, path_counts, reference_data, reference_key, data_by_params):
    # Create widgets
    ladder_selector = widgets.Dropdown(
        options=[(f"Ladder: {s}", s) for s in ladders],
        description='Ladder:',
        style={'description_width': 'initial'}
    )
    
    smoothing_selector = widgets.Dropdown(
        options=[(f"Smoothing: {s}", s) for s in smoothing_params],
        description='Smoothing:',
        style={'description_width': 'initial'}
    )
    
    paths_selector = widgets.Dropdown(
        options=[(f"{p:,} paths", p) for p in path_counts],
        description='Path Count:',
        style={'description_width': 'initial'}
    )
    
    metric_selector = widgets.Dropdown(
        options=[
            ('Price', 'price'), 
            ('S&P Volatility (dPricedVol1)', 'dPricedVol1'),
            ('Apple Volatility (dPricedVol2)', 'dPricedVol2'),
            ('Correlation (dPricedCorr)', 'dPricedCorr'),
            ('S&P Delta (dPricedSpot1)', 'dPricedSpot1'),
            ('Apple Delta (dPricedSpot2)', 'dPricedSpot2'),
            ('S&P Gamma (gamma_spot1)', 'gamma_spot1'),
            ('Apple Gamma (gamma_spot2)', 'gamma_spot2'),
            ('Price Computation Time', 'price_time'),
            ('Risk Computation Time', 'risk_time')
        ],
        description='Metric:',
        style={'description_width': 'initial'}
    )
    
    method_selector = widgets.SelectMultiple(
        options=[('Base', 'base'), ('Smoothed', 'smooth'), ('AADC', 'aadc')],
        value=['base', 'smooth', 'aadc'],
        description='Methods:',
        rows=3,
        style={'description_width': 'initial'}
    )
    
    plot_type = widgets.RadioButtons(
        options=[f'Values vs Ladder'], #, 'Convergence vs Path Count', 'Computation Time', 'RMSE vs Path Count'],
        description='Plot Type:',
        style={'description_width': 'initial'},
        layout={'width': 'max-content'}
    )
    
    #update_button = widgets.Button(description='Update Plot')
    
    # Layout
    controls = widgets.VBox([
        widgets.HBox([ladder_selector, smoothing_selector, paths_selector]),
        widgets.HBox([metric_selector, method_selector]),
        plot_type#,
        #update_button
    ])
    
    output = widgets.Output()
    
    # Function to plot values vs volatility
    def plot_vs_ladder(ladder, smoothing, paths, metric, methods):
        if (smoothing, paths, ladder) not in data_by_params:
            return "Data not available for these parameters."
        
        df = data_by_params[(smoothing, paths, ladder)]
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Add series for selected methods
        if 'base' in methods:
            ax.plot(df[ladder], df[metric], 'o-', linewidth=2, label=f'Base')
        if 'smooth' in methods and f'smooth_{metric}' in df.columns:
            ax.plot(df[ladder], df[f'smooth_{metric}'], 's-', linewidth=2, label=f'Smoothed')
        if 'aadc' in methods and f'aadc_{metric}' in df.columns:
            ax.plot(df[ladder], df[f'aadc_{metric}'], 'x--', linewidth=2, label=f'AADC')
        
        # Reference data
        if metric != 'price_time' and metric != 'risk_time':
            ax.plot(df[ladder], reference_data[ladder][metric], '--', linewidth=2, label=f'BASE Reference ({reference_key[ladder][1]:,} paths)')
#            ax.axhline(y=reference_data[metric].iloc[0], color='r', linestyle='--', 
#                      label=f'Reference ({reference_key[1]:,} paths)')
        
        ax.set_xlabel(ladder)
        ax.set_ylabel(metric_selector.options[metric_selector.index][0])
        ax.set_title(f'{metric_selector.options[metric_selector.index][0]} vs {ladder}\n'
                    f'Smoothing={smoothing}, Paths={paths:,}')
        ax.legend()
        ax.grid(True)
        
        return fig
    
    # Function to plot convergence vs path count
    def plot_convergence(smoothing, metric, methods):
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Reference value for this metric
        ref_value = reference_data[metric].iloc[0]
        
        # Plot data for each path count
        path_data = {}
        for p in path_counts:
            key = (smoothing, p)
            if key in data_by_params:
                path_data[p] = data_by_params[key]
        
        # Sort by path count
        sorted_paths = sorted(path_data.keys())
        
        # Calculate errors against reference
        base_errors = []
        smooth_errors = []
        aadc_errors = []
        
        for p in sorted_paths:
            df = path_data[p]
            
            if 'base' in methods:
                base_error = abs(df[metric].iloc[0] - ref_value) / (abs(ref_value) if ref_value != 0 else 1)
                base_errors.append(base_error)
                
            if 'smooth' in methods and f'smooth_{metric}' in df.columns:
                smooth_error = abs(df[f'smooth_{metric}'].iloc[0] - ref_value) / (abs(ref_value) if ref_value != 0 else 1)
                smooth_errors.append(smooth_error)
                
            if 'aadc' in methods and f'aadc_{metric}' in df.columns:
                aadc_error = abs(df[f'aadc_{metric}'].iloc[0] - ref_value) / (abs(ref_value) if ref_value != 0 else 1)
                aadc_errors.append(aadc_error)
        
        # Plot errors
        if 'base' in methods and base_errors:
            ax.plot(sorted_paths, base_errors, 'o-', linewidth=2, label='Base')
        if 'smooth' in methods and smooth_errors:
            ax.plot(sorted_paths, smooth_errors, 's-', linewidth=2, label='Smoothed')
        if 'aadc' in methods and aadc_errors:
            ax.plot(sorted_paths, aadc_errors, 'x--', linewidth=2, label='AADC')
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Number of Monte Carlo Paths')
        ax.set_ylabel('Relative Error')
        ax.set_title(f'Convergence of {metric_selector.options[metric_selector.index][0]}\n'
                    f'Smoothing={smoothing}, Reference={reference_key[1]:,} paths')
        ax.legend()
        ax.grid(True)
        
        return fig
    
    # Function to plot computation time
    def plot_computation_time(smoothing, time_type, methods):
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Determine which time column to use
        if time_type == 'price':
            base_col = 'price_time'
            smooth_col = 'smooth_price_time'
            aadc_col = 'aadc_price_time'
            title_prefix = 'Pricing'
        else:
            base_col = 'risk_time'
            smooth_col = 'smooth_risk_time'
            aadc_col = 'aadc_risk_time'
            title_prefix = 'Risk Calculation'
        
        # Plot data for each path count
        path_data = {}
        for p in path_counts:
            key = (smoothing, p)
            if key in data_by_params:
                path_data[p] = data_by_params[key]
        
        # Sort by path count
        sorted_paths = sorted(path_data.keys())
        
        # Extract times
        base_times = []
        smooth_times = []
        aadc_times = []
        
        for p in sorted_paths:
            df = path_data[p]
            
            if 'base' in methods and base_col in df.columns:
                base_times.append(df[base_col].iloc[0])
                
            if 'smooth' in methods and smooth_col in df.columns:
                smooth_times.append(df[smooth_col].iloc[0])
                
            if 'aadc' in methods and aadc_col in df.columns:
                aadc_times.append(df[aadc_col].iloc[0])
        
        # Plot times
        if 'base' in methods and base_times:
            ax.plot(sorted_paths, base_times, 'o-', linewidth=2, label='Base')
        if 'smooth' in methods and smooth_times:
            ax.plot(sorted_paths, smooth_times, 's-', linewidth=2, label='Smoothed')
        if 'aadc' in methods and aadc_times:
            ax.plot(sorted_paths, aadc_times, 'x--', linewidth=2, label='AADC')
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Number of Monte Carlo Paths')
        ax.set_ylabel('Computation Time (s)')
        ax.set_title(f'{title_prefix} Time vs Path Count\nSmoothing={smoothing}')
        ax.legend()
        ax.grid(True)
        
        return fig
    
    # Function to plot RMSE across all volatilities for a given path count
    def plot_rmse_vs_paths(smoothing, metric, methods):
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Plot data for each path count
        path_data = {}
        for p in path_counts:
            key = (smoothing, p)
            if key in data_by_params:
                path_data[p] = data_by_params[key]
        
        # Sort by path count
        sorted_paths = sorted(path_data.keys())
        
        # Calculate RMSE against reference for all volatilities
        base_rmse = []
        smooth_rmse = []
        aadc_rmse = []
        
        for p in sorted_paths:
            df = path_data[p]
            
            if 'base' in methods:
                errors = [(x - y)**2 for x, y in zip(df[metric], reference_data[metric])]
                rmse = np.sqrt(np.mean(errors))
                base_rmse.append(rmse)
                
            if 'smooth' in methods and f'smooth_{metric}' in df.columns:
                errors = [(x - y)**2 for x, y in zip(df[f'smooth_{metric}'], reference_data[metric])]
                rmse = np.sqrt(np.mean(errors))
                smooth_rmse.append(rmse)
                
            if 'aadc' in methods and f'aadc_{metric}' in df.columns:
                errors = [(x - y)**2 for x, y in zip(df[f'aadc_{metric}'], reference_data[metric])]
                rmse = np.sqrt(np.mean(errors))
                aadc_rmse.append(rmse)
        
        # Plot RMSE
        if 'base' in methods and base_rmse:
            ax.plot(sorted_paths, base_rmse, 'o-', linewidth=2, label='Base')
        if 'smooth' in methods and smooth_rmse:
            ax.plot(sorted_paths, smooth_rmse, 's-', linewidth=2, label='Smoothed')
        if 'aadc' in methods and aadc_rmse:
            ax.plot(sorted_paths, aadc_rmse, 'x--', linewidth=2, label='AADC')
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Number of Monte Carlo Paths')
        ax.set_ylabel('RMSE (Root Mean Square Error)')
        ax.set_title(f'RMSE of {metric_selector.options[metric_selector.index][0]} vs Path Count\n'
                    f'Smoothing={smoothing}, Reference={reference_key[1]:,} paths')
        ax.legend()
        ax.grid(True)
        
        return fig
    
    # Function to update the plot
    def update_plot(_):
        with output:
            clear_output(wait=True)

            ladder = ladder_selector.value
            smoothing = smoothing_selector.value
            paths = paths_selector.value
            metric = metric_selector.value
            methods = method_selector.value
            plot_choice = plot_type.value
            
            if plot_choice == f'Values vs Ladder':
                result = plot_vs_ladder(ladder, smoothing, paths, metric, methods)
            elif plot_choice == 'Convergence vs Path Count':
                result = plot_convergence(smoothing, metric, methods)
            elif plot_choice == 'Computation Time':
                result = plot_computation_time(smoothing, metric, methods)
            elif plot_choice == 'RMSE vs Path Count':
                result = plot_rmse_vs_paths(smoothing, metric, methods)
            
            if isinstance(result, str):
                print(result)
            else:
                plt.show()
    
    # Connect button to update function
    # update_button.on_click(update_plot)
    # Connect widget changes directly to the update function
    ladder_selector.observe(update_plot, names='value')
    smoothing_selector.observe(update_plot, names='value')
    paths_selector.observe(update_plot, names='value')
    metric_selector.observe(update_plot, names='value')
    method_selector.observe(update_plot, names='value')
    plot_type.observe(update_plot, names='value')    
    
    # Display the widgets
    display(controls, output)
    
    # Initial update
    update_plot(None)
