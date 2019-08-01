# MeerKAT-Cookbook
Recipes for MeerKAT data interaction and processing


## Data archive
Older MeerKAT files are still HDF5 format, but the data rates from radio telescopes are now so large that combining data into single files cause sizes of Giga- to Tera bytes.
These files are to big for standard io-operations.     
Newer data is stored in a more flexible format, and data is accessed and processed as needed.    
New MeerKAT data format is called MeerKAT Visibility Format (MVF)

Cheat sheets and quick introduction usage examples for
[archive](https://github.com/ska-sa/MeerKAT-Cookbook/tree/master/archive)


## KATDAL
Interacting with any of the MeerKAT observation files is made easy by the `katdal` python library   
[katdal repository](https://github.com/ska-sa/katdal)

```
pip search katdal
katdal (0.13)  - Karoo Array Telescope data access library for interacting with data sets in the MeerKAT Visibility Format (MVF)
```

Introductory Jupyter notebooks illustrating some example data interaction and inspection using 
[katdal](https://github.com/ska-sa/MeerKAT-Cookbook/tree/master/katdal)


## CASA

 -fin-
