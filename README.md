# MGCAMB v3.0
## Modified Growth with CAMB
This is the official repository for the MGCAMB v3.0 patch.  Below there are an introduction to the code and the instructions to install 

## Table of contents
*[Introduction](#intro)
    *[Strucutre of the code]
    *[Referencing MGCAMB](#referencing)
    
*[How to install](#how-to-install)
*[How to run](#how-to-run)
    *[Run the code](#run-code)
    *[Run the test suite](#run-tests)
*[What's new](#whats-new)
*[Examples](#examples)
*[Authors List](#author-list)


## Introduction
Modified Growth with CAMB (MGCAMB) is a patch for the Einstein Boltzmann solver CAMB that intrdouces phenomenological modifications of growth along with dynamical Dark Energy. It includes several phenomenological parametrization 


### Structure of the code
The new MGCAMB patch is structured as in the figure.

<p align="center">
<img src="img/MGCAMB_flowchart.png" width="500" title="Distribution of weights and biases updates" />
</p>

### Referencing MGCAMB
If you use MGCAMB for your scientific work, please cite the following papers.

## How to install
To install MGCAMB in your machine simply run
```bash
git clone https://github.com/sfu-cosmo/MGCAMB.git
cd MGCAMB
make camb
```

## How to run

### Run the code
To run MGCAMB, first modify the  ``` params_MG.ini ``` file. Then run


### Run the test suite
If you want to run the test suite, to produce the consistency plots in our paper


## What's new

Alex Zucca: azucca@sfu.ca
