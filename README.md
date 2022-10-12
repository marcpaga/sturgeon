# sturgeon

## Installation

```
python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
pip3 install .
```

## Usage

This program can be used in two modes: `predict` and `watch`.

Once launched, a log file will be created in the `logs` folder. The location of that folder will be dependent from the current working directory (from where the tool is executed).

### Predict

This mode predicts a set of samples in bed file format given a set of models.
Each sample is treated independently.

Example usage with the demo data:
```
# if the environment has not been activated
source venv/bin/activate
sturgeon predict -i demo/ -o demo/results/ -m models/diagnostics.zip -p
```

Parameters:

- `-i`: path to a single bed file or a folder with bed files.
- `-o`: path where the results will be saved.
- `-m`: model file to be used to predict. It can also be several model files.
- `-p`: it will also make barplots of the results (optional).


### Watch

Not Implemeted

## Data formats

### Input

Input files should be in the bed file format. 
Examples of input files can be found in the `demo` folder.

Only the `methylation_call` and `probe_id` columns are obligatory.
In the `methylation_call` column, methylated probes should be indicated with a value of 1, and non-methylated probes with a value of 0. Non-measured probes should not be included.

### Output

Output files are in the csv format.

They contain a column named `number_probes` that indicated the number of measured probes of that sample. The rest of the columns indicate the possible classes in the model. Values indicate the score that the model gave to each class. Higher scores indicate higher confidence in the prediction. Scores cannot be directly linked to likelihood of correctness of the prediction.

## Software compatibility

This tool has only been tested on: 

- `python 3.7` it might work on other versions.