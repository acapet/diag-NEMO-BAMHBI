# diag-NEMO-BAMHBI
Complete NEMO-BAMHBI outputs with a catalogue of diagnostics

To install, enter a terminal an type: 

    git clone git@github.com:acapet/diag-NEMO-BAMHBI.git
    cd diag-NEMO-BAMHBI
    conda env create -f environment.yml

This may take a little while.
Note that you will need to load a conda module from your cluster environment beforehand, which depends on your system.

To use

    conda activate BAMHBI-DIAG
    python diag.py --dir <DIRECTORY>

`<DIRECTORY>` should be a directory including outputs of NEMO-PISCES simulations.

The script will search for any `*ptrc*.nc` files in the directory and issue a new file `*diag*.nc` with the same name structure, nc attributes, and containing the newly computed 2D diagnsotics.

Note that the file and variable requirements depends on the list of requested diagnostics.

The list of diagnostics can be given as an argument, for isntance : 

    python diag.py --dir ./ --diaglist TPPI zchlmax nitracline
    
Each element should be an entry of the diagnostic catalog (see description in the wiki).


