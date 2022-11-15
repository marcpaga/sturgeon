# sturgeon

Sturgeon is a CNS neural network classifier based on the reference dataset published by [Capper et al., 2018](https://doi.org/10.1038/nature26000).

## Installation

From source:
```
python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
pip3 install . --no-cache-dir
```

If you use a pre-compiled binary then there's no necessity for installation.

## Usage

This program has four main utilities:

- `bamtobed`: convert bam files, output of Guppy, to bed files suitable to be used for prediction using this tool.
- `predict`: predict the CNS type from a bed file.
- `livebam`: watch over an output folder where bam files are written live. Then convert and predict them as they come. This is meant to be used during live sequencing.
- `models`: list, add and delete models. Not strictly necessary, as models can be passed to the previous utilities by path.

Once launched, a log file will be created in the `logs` folder. The location of that folder will be dependent from the current working directory (from where the tool is executed).

CRITICAL NOTES:
- Guppy can output bam index files (.bai), do NOT use those, they are not compatible. Either provide them via `samtools` or do not provide them. This program comes with pysam and will create the index files itself.
- Prediction relies on proper mapping of the methylation calls. This program expects that the reads have aligned to the [T2T reference genome v2.0](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz). Alignment to other reference genomes will very likely give incorrect results.


### `bamtobed`

Convert a bam that contains methylation calls into the adequate format (.bed) so that it can be used for prediction.

To obtain bam files from Guppy you can add the following arguments:
```
.guppy_basecaller \
--config dna_r9.4.1_e8.1_modbases_5mc_cg_hac.cfg \
--bam_out \
--align_ref chm13v2.0.fa
```

Do NOT pass the `--index` flag to get index files, these index files are not compatible.
Also remember to pass the T2T reference genome v2.0 as reference for alignment.

Convert bam to bed example with demo data:

```
source venv/bin/activate # <- Not necessary from compiled binaries

sturgeon bamtobed \
-i demo/bam \
-o demo/bam/out
```

In `demo/bam/out` you should have some output files, the most important one is `merged_probes_methyl_calls.bed` as it is the input for the prediction tool. This input contains the aggregation of all methylation calls across ALL bam files. If each bam file is a different sample, then process them independently in different folders.

### `predict`

This mode predicts a set of samples in bed file format given a set of models.
Each bed file is considered to be a single sample and is treated independently.

Example usage with the demo data:
```
source venv/bin/activate # <- Not necessary from compiled binaries

sturgeon predict \
-i demo/ \
-o demo/results/ \
--model-files models/diagnostics.zip \
--plot-results
```

In `demo/results` there should be a `.csv` file for each sample with the scores for each CNS class. There should also be a `.pdf` with a barplot representing the predicted scores.

Values indicate the score that the model gave to each class. Higher scores indicate higher confidence in the prediction. Scores cannot be directly linked to likelihood of correctness of the prediction.

### `livebam`

This program is meant to be used during live basecalling. It watched over a folder and waits for bam files (output of Guppy e.g.) to be written there. Then it processes them as they come. This program expects that all bam files in that folder come from the same sample, therefore the amount of sequencing for that sample increases over time. In this line, each bam file will not be treated independently, but instead they will be added in a cumulative manner. 

Example usage with demo data:
```
source venv/bin/activate # <- Not necessary from compiled binaries

sturgeon livebam \
-i demo/bam \
-o demo/bam/out_live \
--model-files models/diagnostics.zip \
--plot-results

```

The tool needs to be stopped manually because it will wait infinitely for new bam files. In most systems CTRL+C should interrupt and exit the program.

In the output folder there will be a bunch of intermediate files, the most important ones are:

- `predictions_modelname.csv`: which contains the predicted scores for each CNS class. Each row contains the cumulative predicitions, so row 0 are just the predictions for the first bam file, row 1 are the predictions for the first and second bam files combined, ect.
- `predictions_n_modelname.pdf`: these contains barplots for each of the rows in the previous described csv file.
- `predictions_overtime_modelname.pdf`: this contains a plot that describes the change in scores over time. Only classes with an average score >0.1 over time are plotted.
