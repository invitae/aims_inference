# aims_inference

`aims_inference` is a commandline tool for estimating an individual's genetic ancestry from a set of ancestry-informative markers (AIMs). The package comes with a pre-trained classifier that takes a specific set of 45 AIMs as input and outputs a genetic ancestry prediction at the level of continental superpopulation. New classfiers can also be trained using the 1000 genomes project for different sets of AIMs. 

## Installation

The first step in using the tool is to create and activate a local python environment with the required python packages. This can be done several ways (ex. `conda`) but we demonstrate here with `venv`.

First, locally clone the package into a working directory:

`git clone git@github.com:invitae/aims_inference.git`

`cd aims_inference`

Now, create and activate a new python3 environment then install the required packages:

`python3 -m venv ./ve`

`source ve/bin/activate`

`pip install -r requirements.txt`

## Running a trained model

Once the environment is installed, we can do ancestry inference. The required sample input file is in CSV format. The first column, `sample`, gives the sample name. If there is a self-reported ethnicity associated with the sample, it should be included in the sample name like: `[sample-name]_[ethnicity]`. The subsequent columns are labeled with SNP names (ex. `rs3823159`). These should be the same SNPs that the model was trained on. The pre-trained model SNPs are given in `data/kidd_et_al_aims.txt` and were derived from [this paper](https://pubmed.ncbi.nlm.nih.gov/24508742/) by the Kidd lab at Yale. Each genotype should be encoded by a string of alleles like: `C/C`. 

With the python environment activated, ese the following command to to use the pretrained model on a input file called `test_samples.csv`:

`python infer_ancestry.py -i StratGenomEthnicTestsamples.csv -o predicted_ancestries.csv`

The output file will contain one ancestry prediction for each sample in the `test_samples.csv` file. Ancestries are coded like: 

`EUR`=European, `AFR`=African, `EAS`=East Asian, `SAS`=South Asian, `AMR`=Admixed American

## Generating a confusion matrix of ethnicity and ancestry

If self-reported ethnicities are included and appended to the path, we can produce a confusion matrix of the results, relating self-reported ethnicity to predicted ancestry: 

`python infer_ancestry.py -i StratGenomEthnicTestsamples.csv -o predicted_ancestries.csv -cm`

## Retraining a model on a new set of AIMs


