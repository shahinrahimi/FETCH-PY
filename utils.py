import pandas as pd
import os
import datetime

def log(msg: str):
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[FETCH-PY] {current_time} {msg}")
    
    
def save_results(results, output_folder):
    df = pd.DataFrame(results)
    csv_path = os.path.join(output_folder, "results.csv")
    df.to_csv(csv_path, index=False)
    log(f"Results saved to {csv_path}")