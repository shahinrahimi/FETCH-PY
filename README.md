# FETCH-PY
The FETCH-PY repository contains scripts and tools to automate the FETCH (Fluorescence Expression Through Cell Harvesting) analysis process for FCS (Flow Cytometry Standard) data files.

## Summary
This repository aims to streamline the analysis workflow, which includes quality data checks, applying various gating procedures, and calculating the FETCH score. The automation supports contour-based gating, standard deviation gating, and quadrant gating. The final FETCH score is saved to a CSV file, and potential issues with data integrity are flagged. The repository also includes options for generating PDF outputs of the gating steps to assist in debugging and result verification.

## Installation
clone repository
```
git clone https://github.com/shahinrahimi/FETCH-PY.git
cd FETCH-PY
```
create environment
```
conda env create -n fetchenv --file fetchenv.yml
```
activate environment
```
conda activate fetchenv
```

## Usage
To run the pipeline, use the following command:
```
python main.py -f <target_folder> -e <skip_files> -w <overwrite>

```
 - `-f` or `--folder`: Path to the folder containing `.fcs` files for analysis (default is `example`).
 - `-e` or `--skip-files`: List of `.fcs` files to skip (optional).
 - `-w` or `--overwrite`: Whether to overwrite existing results (default is `True`).

 ### Example
```
python main.py -f example -e skip1.fcs skip2.fcs -w True
```

 ## Parameters
 ### Biexponential Transformation Parameters
 The pipeline uses a biexponential transformation for data processing. You can adjust the following parameters in the Config class:
  ```python
 class Config:
    max_value: int = 10000000
    width: int = -100
    pos: float = 4.9
    neg: float = 0

 ```
 - `max_value`: The maximum value for the data to be transformed. Default is `10,000,000`.
 - `width`: The width of the linear region of the transformation. Controls the smoothness of the transition between linear and logarithmic regions. Default is `-100`.
 - `pos`: The positive minimum value for the transformation. Defines where the transition to the linear region starts. Default is `4.9`.
 - `neg`: The negative minimum value for the transformation. Sets the minimum threshold for negative values. Default is `0`.

