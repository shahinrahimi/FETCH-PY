import pandas as pd
import os
def save_fetch_score(fetch_scores, output_folder):
    df = pd.DataFrame(fetch_scores)
    csv_path = os.path.join(output_folder, "fetch_scores.csv")
    df.to_csv(csv_path, index=False)
    print(f"Fetch scores saved to {csv_path}")