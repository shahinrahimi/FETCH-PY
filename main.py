import argparse
import os
import shutil
import flowkit as fk
from gating import check_channels, first_gating_plot, second_gating_plot, third_gating_plot
from utils import save_results, log

def process_files(target_folder, skip_files, overwrite):
    # Check if target folder exists
    if not os.path.exists(target_folder):
        log(f"Target folder {target_folder} does not exist.")
        return
    # define output folder
    output_base_folder = f"{target_folder}_results"
    if not os.path.exists(output_base_folder):
        os.makedirs(output_base_folder)
        log(f"Created results folder {output_base_folder}")
        
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
        
        if check_channels(sample):
            # Create the output folder for the file has required channels
            output_folder = os.path.join(output_base_folder, os.path.splitext(fcs_file)[0])
        
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
            fetch_score = third_gating_plot(df, output_folder)
            results.append({
                'file_name': fcs_file,
                'has_required_channels': True,
                'fetch_score': fetch_score
            })
        else:
            log(f"{fcs_file} is missing required channels, skipping analysis.")
            results.append({
                'file_name': fcs_file,
                'has_required_channels': False,
                'fetch_score': 'N/A'
            })
    
    # Save results
    save_results(results, output_base_folder, f"{target_folder}_results.csv")

    
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