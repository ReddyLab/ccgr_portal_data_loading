# CCGR Portal Data Generation

This repository contains the scripts to generate portal data from public sources. The scripts need to

- Generate experiment and analysis metadata
- Download the data itself (potentially)
- Modify the data to match the portal format

Portal admins can then use a [script included with the portal](https://github.com/ReddyLab/cegs-portal/blob/main/scripts/data_loading/bulk_upload.py) to upload the data to their portal.
