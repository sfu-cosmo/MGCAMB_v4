# MGCAMB v3.0
## Modified Growth with CAMB
This is the official repository for the MGCAMB v3.0 patch.  Below there are an introduction to the code and the instructions to install 

## Table of contents
*[Introduction](#introduction)
   *[Structure of the code](#structure-of-the-code)
   *[Consistency of the code]
   *[Referencing MGCAMB](#referencing)
*[How to install](#how-to-install)
*[How to run](#how-to-run)
   *[Run the code](#run-the-code)
   *[Run the test suite](#run-the-test-tests)
*[What's new](#whats-new)
*[Examples](#examples)
*[Authors List](#authors-list)


## Introduction
Modified Growth with CAMB (MGCAMB) is a patch for the Einstein Boltzmann solver CAMB that intrdouces phenomenological modifications of growth along with dynamical Dark Energy. It includes several phenomenological parametrization 


### Structure of the code
The new MGCAMB patch is structured as in the figure.

<p align="center">
<img src="img/MGCAMB_flowchart.png" width="1000" title="MGCAMB code structure" />
</p>


### Consistency of the code
The General Relativity (GR) limit of the code has been tested 

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

```bash
./camb params_MG.ini
```


### Run the test suite
If you want to run the test suite to produce the consistency plots in our paper, then run

```bash
cd mgcamb_tests

```


## What's new

## Examples

## Authors List
Main Developer:
- Alex Zucca 

Original Code Developer:
-Gong-Bo Zhao
-Alessandra Silvestri
-Levon Pogosian

Alex Zucca: azucca@sfu.ca
