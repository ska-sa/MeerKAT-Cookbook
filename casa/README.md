## Standard MeerKAT L-band flagging
Refer to notebook [MeerKAT-Cookbook/utils/L-band RFI frequency
flagging.ipynb](https://github.com/ska-sa/MeerKAT-Cookbook/blob/master/utils/L-band%20RFI%20frequency%20flagging.ipynb)

Apply standard L-band RFI flags for MeerKAT to CASA measurement set using utility script:
* `casa --log2term --noglogger`
* `run flagging_mkat_lband.py -h`

Usage example:    
```
run flagging_mkat_lband.py --msfile <filename.ms> --zeros --bp --lband --gps --glonass --galileo
--afristar --iridium --inmarsat`
```


## Standard MeerKAT antenna based calibration
* `casa --log2term --noglogger`
* `run calibrate_mkat_lband.py -h`

Usage example:    
```
run calibrating_mkat_lband.py --msfile <filename.ms> -f <f_cal> [-t <target>] [-g <g_cals>] [-b <b_cals>] [-r <ref_ant>]
```

