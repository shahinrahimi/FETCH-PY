import numpy as np
import flowkit as fk
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from scipy import stats


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
    

def first_gating_plot(df) -> plt:
    x_label = Chan.FSC_A
    y_label = Chan.SSC_A 
    x = df[x_label]
    y = df[y_label]
    plt.figure(figsize=(10,10))
    plt.scatter(x, y,alpha=0.2, c='black', s=1)
    sns.kdeplot(data=df, x = x_label, y = y_label,fill=False, alpha=1, levels=25, cmap="coolwarm",linewidths=1)
    # Set axis limits
    plt.xlim(0, 4.2e6)
    plt.ylim(0, 4.2e6)
    plt.title(f'{x_label} vs {y_label} Scatter Plot')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    return plt

def second_gating_plot(df, sd_df=2) -> plt:
    x_label = Chan.FSC_A
    y_label = Chan.FSC_H
    x = df[x_label]
    y = df[y_label]
    # Perform linear regression to get the line of best fit
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    # Calculate the fitted values (the line of best fit)
    fitted_values = slope * x + intercept
    # Calculate the residuals and standard deviation of residuals
    residuals = y - fitted_values
    std_dev = np.std(residuals)
    # Define the upper and lower gate lines (4 standard deviations away)
    upper_gate = fitted_values + sd_df * std_dev
    lower_gate = fitted_values - sd_df * std_dev
    
    # Plot the data and the gating lines
    plt.figure(figsize=(10, 10))
    
    # Scatter plot
    plt.scatter(x, y, alpha=0.5, c='yellow', s=1)
    
    # Plot the line of best fit
    plt.plot(x, fitted_values, color='blue', label='Line of Best Fit')
    
    # Plot the upper and lower gating lines
    plt.plot(x, upper_gate, color='red', linestyle='--', label=f'Upper Gate ({sd_df} SD)')
    plt.plot(x, lower_gate, color='green', linestyle='--', label=f'Lower Gate ({sd_df} SD)')
    
    # Set axis limits if needed
    plt.xlim(0, 4.2e6)
    plt.ylim(0, 4.2e6)
    
    # Labeling the ticks with "5 positive decades"
    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{int(val):e}'))
    plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{int(val):e}'))
    
    # Add labels and legend
    plt.title(f'{x_label} vs {y_label} with Gating Lines')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    
    return plt

def third_gating_plot(sample: fk.Sample) -> plt:
    x_label = Chan.Comp_mEmerald_A
    y_label = Chan.Comp_mCherry_A
    
    biex_xform = fk.transforms.WSPBiexTransform(
        'biex',
        max_value=10000000,
        positive=4.9,
        width=-100,
        negative=0
    )
    
    sample.apply_transform(biex_xform)
    df = sample.as_dataframe(source="xform")
    # Plot the data and the gating lines
    plt.figure(figsize=(10, 10))
    x = df[x_label]
    y = df[y_label]
    # Determine the midpoint for gating (optional: use specific thresholds if needed)
    x_median = np.median(x)
    y_median = np.median(y)
    
    plt.scatter(x,y, alpha=0.5, c='black', s=1)
    plt.axvline(x=x_median, color='black', linestyle='--', linewidth=1)
    plt.axhline(y=y_median, color='black', linestyle='--', linewidth=1)
    
    plt.title(f'{x_label} vs {y_label} Quadrant Gating')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    return plt