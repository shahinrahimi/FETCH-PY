import os
import numpy as np
import flowkit as fk
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from config import Config

plt.switch_backend('agg')
plt.ioff()

class Chan():
    SSC_H = 'SSC-H'
    SSC_A = 'SSC-A'
    FSC_H = 'FSC-H'
    FSC_A = 'FSC-A'
    SSC_B_H = 'SSC-B-H'
    SSC_B_A = 'SSC-B-A'
    Comp_mEmerald_A = 'mEmerald-A'
    Comp_mCherry_A = 'mCherry-A'
    AF_A = 'AF-A'
    Time = 'Time'
    @classmethod
    def required_channels(cls) -> list[str]:
        return [getattr(cls, attr) for attr in dir(cls) if not attr.startswith('__') and not callable(getattr(cls, attr))]
def get_missing_channels(sample: fk.Sample) -> list[str]:
    chns = [ch.strip().upper() for ch in sample.channels["pnn"]]
    required_chns = [ch.strip().upper() for ch in Chan.required_channels()]

    missing_channels = [ch for ch in required_chns if ch not in chns]
    return missing_channels   
   

def first_gating_plot(df: pd.DataFrame, sample_name: str, output_folder: str, target_folder: str) -> pd.DataFrame:
    x_label = Chan.FSC_A
    y_label = Chan.SSC_A 
    x = df[x_label]
    y = df[y_label]
    kde = gaussian_kde(np.vstack([x, y]), bw_method='scott')
    # Create a grid to evaluate the KDE
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    X, Y = np.meshgrid(np.linspace(x_min, x_max, 100),
                       np.linspace(y_min, y_max, 100))
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = kde(positions).reshape(X.shape)
    
    debris_threshold = 25000
    max_fsc_a, max_ssc_a = np.max(x), np.max(y)
    valid_points = (x> debris_threshold) & (y> debris_threshold) & (x < (max_fsc_a-1000)) & (y < (max_ssc_a-1000))
    
    # Plot 1: Scatter plot (points)
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, alpha=0.4, s=1, c='lightgrey', label='Original Data')
    plt.scatter(x[valid_points], y[valid_points], alpha=0.4, s=1, c='black', label='Gated Data (Gate 1)')
    # plt.title(f'{x_label} vs {y_label} Scatter Plot\n{sample_name}')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.savefig(os.path.join(output_folder,f"{target_folder}_{sample_name}_g1_scatter.pdf"), format="pdf")
    
    # Plot 2: Contour plot (KDE)
    plt.figure(figsize=(8, 6))
    plt.contour(X, Y, Z, cmap='coolwarm', levels=25, linewidths=0.5)
    # plt.title(f'{x_label} vs {y_label} KDE Contour Plot\n{sample_name}')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig(os.path.join(output_folder,f"{target_folder}_{sample_name}_g1_contour.pdf"), format="pdf")

    # Plot 3: Combined plot
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, alpha=0.4, s=1, c='lightgrey', label='Original Data')
    plt.contour(X, Y, Z, cmap='coolwarm', levels=25, linewidths=0.5)
    # plt.title(f'{x_label} vs {y_label} KDE Contour Plot\n{sample_name}')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig(os.path.join(output_folder,f"{target_folder}_{sample_name}_g1_combined.pdf"), format="pdf")
    
    filtered_df = df[valid_points]
    return filtered_df

    

def second_gating_plot(df: pd.DataFrame, sample_name: str, output_folder: str, target_folder: str, sd_df=2) -> pd.DataFrame:
    x_label = Chan.FSC_A
    y_label = Chan.SSC_A
    fsc_a = df[Chan.FSC_A].values
    fsc_h = df[Chan.FSC_H].values
    lr = LinearRegression()
    
    # Reshape data for linear regression
    fsc_a_reshaped = fsc_a.reshape(-1, 1)
    fsc_h_reshaped = fsc_h.reshape(-1, 1)
    
    # Fit linear regression model
    lr.fit(fsc_a_reshaped, fsc_h_reshaped)
    # Predict FSC-H values and calculate norms (distances from the fitted line)
    predicted_fsc_h = lr.predict(fsc_a_reshaped).flatten()
    norms = np.abs(fsc_h - predicted_fsc_h)
    
    # Calculate standard deviation of norms and apply the second gate (4 std deviations)
    std_norm = np.std(norms)
    gate_threshold = sd_df * std_norm

    # Apply the second gate by excluding points outside the threshold
    valid_points = norms < gate_threshold
    filtered_data = df[valid_points]
    
    # Plot the data and the gating lines
    plt.figure(figsize=(8, 8))
    plt.scatter(fsc_a,fsc_h, c='lightgrey', s=1, label="Original Data")
    plt.scatter(fsc_a[valid_points], fsc_h[valid_points], c='black',s=1 , label='Gated Data (Gate 2)')
    # Plot the line of best fit
    plt.plot(fsc_a, predicted_fsc_h, color='green', label='Fitted Line', linewidth=2)
    upper_bound = predicted_fsc_h + gate_threshold
    lower_bound = predicted_fsc_h - gate_threshold
    plt.plot(fsc_a, upper_bound, color='red', linestyle='--', label='Upper Bound (+4σ)', linewidth=1)
    plt.plot(fsc_a, lower_bound, color='red', linestyle='--', label='Lower Bound (-4σ)', linewidth=1)
    
    # Labeling the ticks with "5 positive decades"
    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{int(val):e}'))
    plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{int(val):e}'))
    
    plt.xticks(rotation=45)
    
    # Add labels and legend
    # plt.title(f'{x_label} vs {y_label} with Gating Lines\n{sample_name}')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    
    plt.savefig(os.path.join(output_folder,f"{target_folder}_{sample_name}_g2.pdf"), format="pdf")
    
    return filtered_data

def third_gating_plot(df: pd.DataFrame, sample_name: str ,output_folder: str, target_folder: str) -> float | str | None:
    x_label = Chan.Comp_mEmerald_A
    y_label = Chan.Comp_mCherry_A

    sam = fk.Sample(df, sample_id="Gated Test")
    df2 = sam.as_dataframe(source='raw')
    x = df2[Chan.Comp_mEmerald_A]
    y = df2[Chan.Comp_mCherry_A]
    
    # specify transform values
    max_value = max(x.max(), y.max())
    rng = max(x.max() - x.min(), y.max() - y.min())
    width = rng * Config.width_multi
    biex_xform = fk.transforms.WSPBiexTransform(
        'biex',
        max_value=max_value,
        positive=Config.pos,
        width=-width,
        negative=Config.neg
    )
    
    sam.apply_transform(biex_xform)
    df2 = sam.as_dataframe(source='xform')
    
    # new transformed value
    x = df2[x_label]
    y = df2[y_label]
    # old untransformed values
    xx = df[x_label]
    yy = df[y_label]
    
    # Define autofluorescence cutoff
    autofluorescence_cutoff = 4200
    
    # Kernel Density Estimation
    data = np.vstack([x, y]).T
    # Grid search over bandwidth in log-space
    bandwidths = np.logspace(-2, 1, 50)
    grid = GridSearchCV(KernelDensity(), {'bandwidth': bandwidths}, cv=5)
    grid.fit(data)
    best_bandwidth = grid.best_estimator_.bandwidth_

    kde = KernelDensity(bandwidth=best_bandwidth).fit(data)

    # Filter density scores based on autofluorescence cutoff
    valid_data = (x < autofluorescence_cutoff) & (y < autofluorescence_cutoff)
    kde_values = kde.score_samples(np.vstack([x[valid_data], y[valid_data]]).T)
    
    # Get the contour levels
    number_of_levels = 30
    levels = np.linspace(0, kde_values.max(), number_of_levels)
    # normalize levels [0,1]
    min_val = levels.min()
    max_val = levels.max()
    
    normalize_levels = (levels - min_val) / (max_val - min_val)
    sorted_normalize_levels = np.sort(normalize_levels)
    
    plt.figure(figsize=(10, 10))
    ax = sns.kdeplot(
        x=x, 
        y=y,
        levels=sorted_normalize_levels, 
        bw_method="scott",
        cmap="coolwarm", 
        fill=True, 
        thresh=0.001, 
        linewidths=0.1,
        legend=True
    )
    
    def get_first_best_lines(paths, lvl):
        for j, path in enumerate(paths):
            level = sorted_normalize_levels[j % len(sorted_normalize_levels)]
            if level < lvl:
                continue
            vertices = path.vertices
            x_mean = np.mean(x)
            y_mean = np.mean(y)
            x_vertices = vertices[:, 0]
            y_vertices = vertices[:, 1]
            mask = (x_vertices < x_mean) & (y_vertices < y_mean)  # grab Q4 points
            valid_x_vertices = x_vertices[mask]
            valid_y_vertices = y_vertices[mask]
            vline = np.max(valid_x_vertices)
            hline = np.max(valid_y_vertices)
            return vline, hline
        
    target_contour = ax.collections[0]  # collection has one contour level
    paths = target_contour.get_paths()
    vline, hline = get_first_best_lines(paths, 0.6)
    
    # Quadrant calculations
    quadrant_1 = df2[(df2[x_label] < vline) & (df2[y_label] > hline)]  # HIGH Y and LOW X (High mCherry)
    quadrant_3 = df2[(df2[x_label] > vline) & (df2[y_label] < hline)]  # LOW Y and HIGH X (High eEmerald)
    quadrant_2 = df2[(df2[x_label] > vline) & (df2[y_label] > hline)]  # HIGH Y and HIGH X (Double Transfected)
    quadrant_4 = df2[(df2[x_label] < vline) & (df2[y_label] < hline)]  # LOW Y and LOW X (Untransfected)
    
    R = len(quadrant_1)  # mCherry cells (Red)
    G = len(quadrant_3)  # mEmerald cells (Green)
    D = len(quadrant_2)  # Double transfected cells
    U = len(quadrant_4)  # Untransfected cells
    total_population = D + G + R + U
    
    # Calculate fetch score
    fetch_score = D / (D + G + R) if (D + G + R) != 0 else None
    
    
    # Plotting the quadrants and fetch score
    xlim = plt.gca().get_xlim()
    ylim = plt.gca().get_ylim()
    bbox = dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5')

    # Disqualify based on untransfected proportion or fetch score
    if total_population != 0:
        untransfected_proportion = U / total_population
        if untransfected_proportion > 0.9:
            fetch_score = 0
        elif fetch_score is not None and fetch_score > 0.9:
            fetch_score = 0
    # Plot 1: Counter Plot (KDE)
    plt.figure(figsize=(10, 10))
    sns.kdeplot(
        x=x, 
        y=y,
        levels=sorted_normalize_levels, 
        bw_method="scott",
        cmap="coolwarm", 
        fill=False, 
        thresh=0.001, 
        linewidths=1,
        legend=True
    )
    plt.axhline(y=hline, color='black', linestyle='--', linewidth=1)
    plt.axvline(x=vline, color='black', linestyle='--', linewidth=1)
    plt.text(xlim[1], ylim[1], f'Q2 D: {D}', fontsize=10, verticalalignment='top', horizontalalignment='right', color='blue', bbox=bbox)
    plt.text(xlim[0], ylim[1], f'Q1 R: {R}', fontsize=10, verticalalignment='top', horizontalalignment='left', color='red', bbox=bbox)
    plt.text(xlim[1], ylim[0], f'Q3 G: {G}', fontsize=10, verticalalignment='bottom', horizontalalignment='right', color='green', bbox=bbox)
    plt.text(xlim[0], ylim[0], f'Q4 U: {U}', fontsize=10, verticalalignment='bottom', horizontalalignment='left', color='black', bbox=bbox)
    # plt.title(f'{x_label} vs {y_label} Contour Plot\n{sample_name}')
    plt.savefig(os.path.join(output_folder,f'{target_folder}_{sample_name}_g3_contour.png'))
    
    # Plot 2: Scatter Plot (Data points)
    plt.figure(figsize=(10, 10))
    plt.scatter(x, y, alpha=0.1, c='black', s=1, label='Data Points')
    plt.axhline(y=hline, color='black', linestyle='--', linewidth=1)
    plt.axvline(x=vline, color='black', linestyle='--', linewidth=1)
    plt.text(xlim[1], ylim[1], f'Q2 D: {D}', fontsize=10, verticalalignment='top', horizontalalignment='right', color='blue', bbox=bbox)
    plt.text(xlim[0], ylim[1], f'Q1 R: {R}', fontsize=10, verticalalignment='top', horizontalalignment='left', color='red', bbox=bbox)
    plt.text(xlim[1], ylim[0], f'Q3 G: {G}', fontsize=10, verticalalignment='bottom', horizontalalignment='right', color='green', bbox=bbox)
    plt.text(xlim[0], ylim[0], f'Q4 U: {U}', fontsize=10, verticalalignment='bottom', horizontalalignment='left', color='black', bbox=bbox)
    # plt.title(f'{x_label} vs {y_label} Scatter Plot\n{sample_name}')
    plt.savefig(os.path.join(output_folder,f'{target_folder}_{sample_name}_g3_scatter.png'))
    
    # Plot 3: Combined Plot
    plt.figure(figsize=(10, 10))
    sns.kdeplot(
        x=x, 
        y=y,
        levels=sorted_normalize_levels, 
        bw_method="scott",
        cmap="coolwarm", 
        fill=False, 
        thresh=0.001, 
        linewidths=1,
        legend=True
    )
    plt.scatter(x, y, alpha=0.1, c='black', s=1, label='Data Points')
    plt.axhline(y=hline, color='black', linestyle='--', linewidth=1)
    plt.axvline(x=vline, color='black', linestyle='--', linewidth=1)
    plt.text(xlim[1], ylim[1], f'Q2 D: {D}', fontsize=10, verticalalignment='top', horizontalalignment='right', color='blue', bbox=bbox)
    plt.text(xlim[0], ylim[1], f'Q1 R: {R}', fontsize=10, verticalalignment='top', horizontalalignment='left', color='red', bbox=bbox)
    plt.text(xlim[1], ylim[0], f'Q3 G: {G}', fontsize=10, verticalalignment='bottom', horizontalalignment='right', color='green', bbox=bbox)
    plt.text(xlim[0], ylim[0], f'Q4 U: {U}', fontsize=10, verticalalignment='bottom', horizontalalignment='left', color='black', bbox=bbox)
    # plt.title(f'{x_label} vs {y_label} Combined Plot\n{sample_name}')
    plt.savefig(os.path.join(output_folder,f'{target_folder}_{sample_name}_g3_combined.png'))
    
    return fetch_score