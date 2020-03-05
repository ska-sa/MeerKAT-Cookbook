# MeerKAT archive quick info

The obvious first step is to get access to MeerKAT data via the MeerKAT archive, following the steps
summarised on the repository entry page.
It is assumed the user has some familiarity with the archive user manual.

In general, if no search criteria is specified all data available to the user will be listed.
This include publicly available engineering and commissioning observation data, as well as science
proposal observations.    
**Note:** The biggest difference is that engineering and commissioning observation data are generally
publicly available, while access to science observations are under restrictions as documented in the
MeerKAT data access policy.

Important terminology that will pop up:

* The **Observer Proposal ID** is assigned when an observation proposal gets submitted. It becomes the
  primary reference for all observations associated with the observations related to the proposal.    
  The assignment of a `proposal ID` is mandatory, with two major differences: a `proposal ID` for
  submitted science observations are assigned following MeerKAT observation request policies.
  Science observation IDs can be identified by the `SCI` prefix.    
  For engineering and commissioning observations, the `proposal ID` is assigned on an ad-hoc basis.
  These observations do not contain any relevant science data, but are useful for training and
  investigation of MeerKAT data and performance.

* The second import ID tag is the **Schedule Block ID**.
  This ID is automatically assigned by the telescope per observation on execution.
  It distinguish the various observations performed, and allows easy search for newer observations in
  larger projects.

  If neither the `proposal ID` or the relevant `observation ID` is specified, all files available are
  listed.


Notebooks presented in this folder provides examples for data interaction using MeerKAT provided
software `katdal`, as well as conversion of MVF files to CASA measurement sets for standard
processing.


 -fin-
