# CCGR Portal Data Generation, Validation and Loading

This repository contains the scripts to generate portal data from public sources, validate said data, and load it (or any other valid data) into in the portal. The generation scripts need to

- Generate experiment and analysis metadata
- Download the data itself (potentially)
- Modify the data to match the portal format (potentially)

The generated data (or any data) can be validated to work with the portal using the `ccgr_ev` (experiment data validation), `ccgr_av` (analysis data validation), and `ccgr_mv` (experiment and analysis metadata validation) scripts.

Portal admins can then use the `ccgr_ul` command to upload data to their portal.

## Installation

`pip install git+https://github.com/ReddyLab/ccgr_portal_data_loading.git`

## Generating Formatted Data for the Portal

`gen_engreitz`, `gen_encode_myc`, `gen_reddylab` can all be used on specific datasets to extract the necessary
information for portal use, as well as generate the experiment/analysis metadata.

## Leftovers

`ccgr_lo` can be used to lift over datasets from one assembly to another.
