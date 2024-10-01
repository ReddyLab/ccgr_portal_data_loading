# CCGR Portal Data Generation, Validation and Loading

This repository contains the scripts to generate portal data from public sources, validate said data, and load it (or any other valid data) into in the portal. The generation scripts need to

- Generate experiment and analysis metadata
- Download the data itself (potentially)
- Modify the data to match the portal format (potentially)

The generated data (or any data) can be validated to work with the portal using the `ccgr_ev`, `ccgr_av`, and `ccgr_mv` scripts.

Portal admins can then use the `ccgr_ul` command to upload data to their portal.

## Installation

`pip install git+https://github.com/ReddyLab/ccgr_portal_data_loading.git`
