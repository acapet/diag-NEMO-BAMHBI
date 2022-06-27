"""
    The script will list all instances of *ptrc*.nc files and issue additional *diag*nc files with same filename structure and attributes.
    
    Args:
        dir (str): directoy of NMO-PISCES outputs (default : './')
        diaglist (list of str): identifiers of required diagnostics. Each should correspond to an entry in the diagnsotic dictionnary. 
""" 
   # The main arguments is a directory.
# The script will list all instances of *ptrc*.nc files and issue additional *diag*nc files with same filename structure and attributes.


import argparse
from genericpath import isfile
import os
from glob import glob
import xarray as xr
import DiagFunctions_NEMOPISCES as diag

# Arguments management #
parser = argparse.ArgumentParser()
parser.add_argument("-p","--printlist", help="Just print the list of available diagnostic definitions and required variables", action="store_true")
parser.add_argument("-v","--verbose", help="increase output verbosity", action="store_true")
parser.add_argument('-d','--dir', type=str, default='./',
 help='Directory containing the NEMO-PISCES outputs. \n Variable requirements depends on diaglist, but we expect "ptrc","gridT", and "diad" files with the same filename structure')
parser.add_argument('-k','--key', type=str, default='',
 help='Key required to be present in NEMO-PISCES ouptut filenames. If given the script will look for "*ptrc*KEY*nc" or "*KEY*ptrc*nc" files instead of simply "*ptrc*nc" ')
parser.add_argument('-l','--diaglist', nargs='+', default=['TPPI','pocF100','TrophicEfficiencyI','zchlmax','CHLatMAX','ratioLargeM','CHLI','CHLM5','OMZextent','pocM200'],
 help= 'List of the diagnostics to be computed and stored in new *diag*.nc files.(eg. " ... -l TPPI pocF100 ")')

args = parser.parse_args()
indir =args.dir
key =args.key
dlist=args.diaglist

if args.verbose:
    print("verbose: on")

if args.printlist:
    diag.diaglist()
    exit()
    

print('Selected Diags : ')
diag.diaglist(dlist)



# Setting up file lists #
flist_p = glob(indir + '*ptrc*'+key+'*nc')
if not flist_p:
    print("Couldn't any files of the form :"+indir + "*ptrc*"+key+"*nc")
    flist_p = glob(indir +'*' +key+'*ptrc*nc')
    if not flist_p:
            print("Couldn't any files of the form :"+indir +'*' +key+'*ptrc*nc')

flist_d = [ f.replace('ptrc','diad') for f in flist_p ] 
flist_g = [ f.replace('ptrc','gridT') for f in flist_p ] 
flist_o = [ f.replace('ptrc','diag') for f in flist_p ] 


if args.verbose:
    print('Will take care of files : ')
    print(flist_p)

#TODO test for alternate files (gridT, diad), depending on dependencies and provide meaningfull error message

for i, (fp, fd, fg, fo) in enumerate(zip(flist_p, flist_d, flist_g, flist_o) ):
    # 1: 'ptrc' files 
    x_p=xr.load_dataset(fp)
    xl=[x_p]
    # 2: 'diad' files 
    if isfile(fd):
        x_d=xr.load_dataset(fd)
        #FIXME TPP is now provided as time instant.
        # For now I just overwrite and assume time_centered coordinates instead.
        x_d=x_d.assign(time_counter=x_p['time_counter'].data )
        x_d=x_d.assign(time_counter_bounds=(('time_counter', 'axis_nbounds'), x_p['time_counter_bounds'].data) ) #FIXME MANIP on origin of time_counter
        xl.append(x_d)
    # 3: 'gridT' files 
    if isfile(fg):
        x_g=xr.load_dataset(fg)
        xl.append(x_g)

    # Ensure we got all we may need
    x_a=xr.merge(xl)


    # Need to define a cell height variable to use xgcm 
    ### 16032022 AC - Had to replace 'depth_bounds' with x_a['deptht'].attrs['bounds']. 
    ### Maybe similar handling will be needed for other variable names, in which case it should be done in a more organized way
    x_a['h']=(x_a[x_a['deptht'].attrs['bounds']][:,1:]-x_a[x_a['deptht'].attrs['bounds']][:,:-1]).squeeze()
    x_a['h'].attrs={'units'     : 'm',
              'long_name' : 'cells height',
              'valid_min' : -1e20,
              'valid_max' : 1e20,
              'cell_methods' : 'time: mean',
              'coordinates': 'lon lat deptht'}

    bibi=diag.add2D(x_a,dlist, verbose=args.verbose)

    #FIXME Make the following cleaner.
    if 'time_centered' in x_a.keys():
        dlo=dlist+['time_centered']
    else:
        dlo=dlist
        
    bibi[dlist].to_netcdf(fo)

    print(fo + ' completed')
    del x_a