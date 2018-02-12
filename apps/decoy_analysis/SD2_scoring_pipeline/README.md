# Computational pipeline for scoring protein designs using structure- and sequence-based metrics

The goal of this project is to develop a computational pipeline for using `PyRosetta` to score protein designs with a variety of structure- and sequence-based metrics.

The long-term goal is to develop pipelines for scoring both monomeric proteins and protein-protein interfaces. Currently, this project only has a pipeline for scoring monomeric designs. These metrics include many of the `Rosetta`-based metrics used in [Rocklin et al., 2017, Science](http://science.sciencemag.org/content/357/6347/168/tab-figures-data).

## Organization

Currently, there is a single Jupyter notebook called `score_monomeric_designs.ipynb` that describes a computational pipeline for scoring monomeric protein designs.

This project also has a variety of subdirectories described later in this README in more detial:

    * `./data/` : input data for the analysis

    * `./scripts/` : scripts used in the analysis 
    
    * `./results/` : results of the analysis

## Data

* `./data/rocklin_rd4_examples/` : A subset of 10 designs from the Rocklin et al. study used as example input in `score_monomeric_designs.ipynb`. These designs were selected from Round 4 of design in the Rocklin study.

* `./data/Rocklin_rd4_relax_scored_filtered_betanov15.sc` : Scores for the Round 4 designs from Rocklin et al. (`./Rocklin_rd4_PDB_files/`), used in [`../test_score_protocols.ipynb`](test_score_protocols.ipynb) to make sure our scoring functions are working as expected. Specifically, I obtained these values from the paper's [Supplementary Materials](http://science.sciencemag.org/content/suppl/2017/07/12/357.6347.168.DC1?_ga=2.236704449.1555020638.1510690697-778426272.1509735927) by downloading `Files S3`. I then copied the file from the following path: `./aan0693_SI_datasets/design_structural_metrics/rd4_relax_scored_filtered_betanov15.sc`

## Scripts

* `./xmls/` : Contains a series of `RosettaScripts` XML files for scoring designs. There is one file per scoring metric.

## Results

* `./results/scores.csv` : Shows the socres for the set of 10 designs from Rocklin et al., as computed using the pipeline in `score_monomeric_designs.ipynb`.

## To do:

* Develop a computational pipeline for scoring protein-protein interfaces
* Expand the computational pipeline for scoring monomeric complexes to include custom metrics from the Rocklin et al. study.