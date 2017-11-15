[![Build Status](https://travis-ci.org/zooniverse/decals.svg?branch=master)](https://travis-ci.org/zooniverse/decals)
[![codecov](https://codecov.io/gh/zooniverse/decals/branch/master/graph/badge.svg)](https://codecov.io/gh/zooniverse/decals)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)


# Decals

Morphological classifications for DECaLs DR5 by the Galaxy Zoo - DECaLS collaboration.

Written by Mike Walmsley.
Based on [previous work](https://github.com/willettk/decals)
by Kyle Willett and Coleman Krawcyzk for DR1 and DR2 respectively.

#### Key Aims

This repo collects the code related to the DECaLS classification project on Galaxy Zoo.
1. Identifying galaxies imaged in both SDSS and DECaLS
2. Downloading the telescope data
3. Creating RGB images for classification
4. Identifying galaxies not yet uploaded from DECaLS DR1 and DR2
5. Aggregating and interpreting the volunteer classifications

Currently (10.11.2017) aggregation (5) is not yet implemented.

#### Files Required

1. The NASA-Sloan Atlas (NSA) catalog v1_0_1. Available from
SDSS DR13 [here](http://skyserver.sdss.org/dr13/en/help/browser/browser.aspx#&&history=description+nsatlas+U).
2. Brick (image tile) catalogs for DECaLS DR5. Available as two separate tables:
[survey-bricks.fits](http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr5/survey-bricks.fits.gz) (all bricks) and [survey-bricks-dr5.fits](http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr5/survey-bricks-dr5.fits.gz) (imaged bricks).
These tables are combined during setup.
Descriptions are [here](http://legacysurvey.org/dr5/files/).
3. All previous DECaLS Galaxy Zoo subjects and associated metadata.
All previous subjects are available as a .csv data dump on request to the Galaxy Zoo team.

More information is available about DECaLS DR5 [here](http://legacysurvey.org/dr5/).

#### Setup

The code involves quite a few files, and it needs to know where to look
 for or save each one.

To use the default settings, create the following directories.
Download the required files with the links above, and place them
in the directories as shown.

If you only want to download DR5:

+ catalogs (dir)
    - NSA_v1_0_1.fits.
    - survey-bricks.fits
    - survey-bricks-dr5.fits
+ subjects (dir)
    - galaxy_zoo_subjects.csv
+ fits (dir)
    + dr5 (dir)
+ jpeg (dir)
    + dr5 (dir)

Or, to enable downloading DR1, DR2 and DR5:

+ catalogs (dir)
    - NSA_v{e.g. 1_0_1}.fits
    - survey-bricks.fits
    - survey-bricks-dr5.fits
    - decals-bricks.fits
    - decals-bricks-dr1.fits
    - decals-bricks-dr2.fits
+ subjects (dir)
    - galaxy_zoo_subjects.csv
+ fits (dir)
    + dr5 (dir)
    + dr2 (dir)
    + dr1(dir)
+ jpeg (dir)
    + dr5 (dir)
    + dr2 (dir)
    + dr1(dir)

If you prefer to set up your own directories and filenames,
you can configure those settings (and all other settings) in the main file [here](https://github.com/zooniverse/decals/python/get_decals_images_and_catalogs.py).


#### Execution
Navigate to the repo root and run

`pip install -r requirements.txt`

This will install all required dependencies.
If you are using anaconda, you should already have most/all of them.

Open python/get_decals_images_and_catalogs
1. Double-check all the filenames point where you've downloaded files
2. Choose the run configuration, as described inline.

Run unit tests if desired. From repo root directory:

`pytest`

Run everything!

`python/get_decals_images_and_catalogs`


#### Methods

Below is a diagram showing the steps performed. Filenames are the program defaults for a DR5 download.

![flowchart](https://www.lucidchart.com/publicSegments/view/7a2db1ef-4fda-460a-9397-1cea25bd9451/image.jpeg)

#### Known Issues

There are no known issues, but plenty of room to improve performance.
