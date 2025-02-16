import argparse
import os
import shutil
import flowkit as fk
import pandas as pd
from gating import check_channels, first_gating_plot, second_gating_plot, third_gating_plot
from utils import save_results, log

def process_files(target_folder, skip_files, overwrite):
    # Check if target folder exists
    if not os.path.exists(target_folder):
        log(f"Target folder {target_folder} does not exist.")
        return
    # list all .fcs files in target folder
    fcs_files = [f for f in os.listdir(target_folder) if f.endswith('.fcs')]
        
    results = []

    for fcs_file in fcs_files:
        if fcs_file in skip_files:
            log(f"Skipping {fcs_file}")
            continue
        # Load the sample 
        filepath = os.path.join(target_folder, fcs_file)
        sample = fk.Sample(filepath)
        # df_events
        df = sample.as_dataframe(source='raw')
        
        # Log available channel names for verification
        log(f"Available channels for {fcs_file}: {list(df.columns)}")
        
        if check_channels(sample):
            # Create the output folder for the file has required channels
            output_folder = os.path.join(target_folder, os.path.splitext(fcs_file)[0])
        
            if os.path.exists(output_folder):
                if overwrite:
                    shutil.rmtree(output_folder)
                    os.makedirs(output_folder)
                    log(f"Overwriting results in {output_folder}")
                else:
                    log(f"Skipping {fcs_file} as overwrite is disabled and folder exists.")
                    continue
            else:
                os.makedirs(output_folder)
            
            # The gating procedures and save the plots
            # should return df by applying the gate
            df = first_gating_plot(df, output_folder)
            # should return df by applying 2nd gate
            df = second_gating_plot(df, output_folder)
            result = third_gating_plot(df, output_folder)
            
            if isinstance(result, tuple):
                fetch_score, gate_boundaries = result
            else:
                fetch_score = result
                gate_boundaries = None
            
            # Dubious check conditions
            dubious = 'no'
            if fetch_score is None or fetch_score > 0.90:
                dubious = 'yes'
                log(f"{fcs_file} marked as dubious due to high FETCH score or missing data.")
                fetch_score = 0
            
            # Calculate red-to-green ratio (r_g) if applicable
            if 'mCherry-A' in df.columns and 'mEmerald-A' in df.columns:
                red_cells = len(df[df['mCherry-A'] > 4200])
                green_cells = len(df[df['mEmerald-A'] > 4200])
                r_g = red_cells / green_cells if green_cells > 0 else float('nan')
            else:
                log(f"{fcs_file} is missing required fluorescence channels for red-to-green ratio calculation.")
                r_g = float('nan')
            
            # Additional dubious checks
            if r_g >= 2 or r_g <= 0.5 or len(df) < 500 or pd.isna(r_g):
                dubious = 'yes'
                log(f"{fcs_file} marked as dubious due to r_g ratio or low cell count.")
                fetch_score = 0
            
            # Add gate boundaries for verification
            if gate_boundaries:
                log(f"{fcs_file} gate boundaries: Vertical Line - {gate_boundaries['vline']}, Horizontal Line - {gate_boundaries['hline']}")
                if gate_boundaries['hline'] == (df['mCherry-A'].max() + df['mCherry-A'].min()) / 2:
                    log(f"Warning: Horizontal gate boundary for {fcs_file} is set at the midpoint of the y-axis. Please verify gating logic.")
            
            results.append({
                'file_name': fcs_file,
                'has_required_channels': True,
                'fetch_score': fetch_score,
                'dubious': dubious
            })
        else:
            log(f"{fcs_file} is missing required channels, skipping analysis.")
            results.append({
                'file_name': fcs_file,
                'has_required_channels': False,
                'fetch_score': 'N/A',
                'dubious': 'N/A'
            })
    
    # Save results
    save_results(results, target_folder)

    
def main():
    parser = argparse.ArgumentParser(description="FETCH Analysis Pipeline")
    parser.add_argument('-f', '--folder', type=str, default="example", help="Path to generate results (default: 'example').")
    parser.add_argument('-e', '--skip-files',type=str, nargs="*", default=[], help="List of .fcs files to skip (default: empty list).")
    parser.add_argument('-w', '--overwrite',type=bool, default=True, help="Whether to overwrite existing results (default: True).")
    
    args = parser.parse_args()
    target_folder = args.folder
    skip_files = args.skip_files
    overwrite = args.overwrite
    # process the files bases on the provided argument
    process_files(target_folder, skip_files, overwrite)

if __name__ == "__main__":
    main()
