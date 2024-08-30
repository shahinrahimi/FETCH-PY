import argparse
import os
import shutil
from gating import first_gating_plot, second_gating_plot, third_gating_plot
import flowkit as fk

def process_files(target_folder, skip_files, overwrite):
    # Ensure the target folder exists
    if not os.path.exists(target_folder):
        print(f"Target folder {target_folder} does not exist.")
        return
    # list all .fcs files in target folder
    fcs_files = [f for f in os.listdir(target_folder) if f.endswith('.fcs')]

    for fcs_file in fcs_files:
        if fcs_file in skip_files:
            print(f"Skipping {fcs_file}")
            continue
        # Create the output folder for the file
        output_folder = os.path.join(target_folder, os.path.splitext(fcs_file)[0])
        
        if os.path.exists(output_folder):
            if overwrite:
                shutil.rmtree(output_folder)
                os.makedirs(output_folder)
                print(f"Overwriting results in {output_folder}")
            else:
                print(f"Skipping {fcs_file} as overwrite is disabled and folder exists.")
                continue
        else:
            os.makedirs(output_folder)
            
        # Perform the analysis and save result
        filepath = os.path.join(target_folder, fcs_file)
        sample = fk.Sample(filepath)
        df = sample.as_dataframe(source='raw')
        plt = first_gating_plot(df)
        plt.savefig(os.path.join(output_folder, "first_gating.pdf"), format="pdf")
        plt = second_gating_plot(df)
        plt.savefig(os.path.join(output_folder,"second_gating.pdf"), format="pdf")
        plt = third_gating_plot(sample)
        plt.savefig(os.path.join(output_folder,"third_gating.pdf"), format="pdf")
    
    
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