<img src="dep_sign.png" width=120, height=120 align="left" />

# BEvis

BEvis is Bacteria Exchanges Visualization. Here we provide a script to facilitate with the visualisation of source compound processing via bacterial community metabolism as described in Yusim, J et al. 2024.
This script utilizes results from python based metabolic models processing tools, Cobrapy and COMETS to generate a network visualisation.
The script uses networkx package for network visualisation.
The network is constructed using layers. These layers are associated with user-defined group of bacteria and metabolites.
This layered network facilitates with visualisation of source compound processing via bacterial community metabolism.
Visualisation requires associating graphical features, e.g., edge thickness, vertex size, with numerical features, such as flux, biomass.
Therefore, user defined parameters control for the fine-tuning of graphical features.
Modification of user provided parameters is performed within a single config file, BEvis_config.yml. 

This tool accepts input files generated from both KBASE and CARVME (BIGG format).
Therefore the following must confirm to single format (KBASE or BIGG):

1. BEvis_config.yml
2. Compounds as in the format of metabolic model used to generate input files for BEvis.
e.g., "cpd00030_e0" for KBASE and glc__D_c (BIGG format, gereneated using CAREVME)

The provided example uses input files in KBASE

In addition, bacteria names must be consistent within input files (including file names) and  BEvis_config.yml.

## Input files

- Bacteria biomass in time (df_BM.csv).
This is a csv file with "cycle","species","biomass" as header.

- In separate folder, data/model_exchanges, exchange files of each model(e.g, Mesotoga_infera.csv) 
exchanges of each compound in specific cycle.

- A dictionary of: compound id, compound name and formula (dicComp_df.csv)

- Metabolites amounts in time (df_metabolites.csv)

## Determine cycle (day) for visualization

Cycle (day) for visualization, c\_i is calculated as fixed number of cycles , nc, prior to c_max. c_i = c_max - nc.

Where c_max is the cycle which maximal growth (biomass = 0.1) is achieved and nc is user provided.

Lower c_i, or higher nc, corresponds to lower biomass point, during the bacterial exponential growth phase.

## Dependencies

* [matplotlib==3.2.2]
* [networkx==3.1]
* [numpy==1.23.5]
* [pip==22.3.1]
* [python==3.8.0]

### pip:

* [pandas==1.5.1]
* [pyyaml==6.0]

## Installation

### Download and install BEvis on linux (Ubuntu 20.04)

Clone repository using git clone

or downlowd zip and extract repository 

### Create virtual environment and install dependencies

```shell
# Create venv and install dependencies #

conda env create -f BEvis_env.yml

```

### Activate virtual environment and move to working directory  

```shell

conda activate BEvis_env

cd BEvis

```

#### run script

```shell

cd ./data/model_exchanges/

ls *.csv > ../model_id_list.txt 

cd  ../../

python BEvis.py

#Or a single bash script:

bash BEvis_all.sh


```


## Contributors

[Shlomit Medina](https://www.freilich-lab.com/shlomit-medina )
[Gon Carmi](https://www.freilich-lab.com/members) \
[Shiri Freilich](https://www.freilich-lab.com/shiri-detailes )
[Raphy Zarecki](https://www.linkedin.com/in/raphy-zarecki-3412663/?originalSubdomain=il)
[Jenny Yusim](https://www.freilich-lab.com/jenny-details) \

## References

Yusim, J et al. 2024.

## Funding

This work was funded by the Joint NSFC-ISF Research Grant [grant number 3164/19]
