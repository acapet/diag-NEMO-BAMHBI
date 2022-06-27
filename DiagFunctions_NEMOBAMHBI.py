from numpy import append
import xarray as xr
import xgcm

# Catalogue of diagnostics for PISCES outputs
# A. Capet - acapet@uliege.be - Feb 2022

ddiag2D = { 'NPPO'      :  { 'req' : ['NPP','PhytoNitrateReduction'] , 
                            'attrs' : {'units'     : 'mmol O2 m-3 s-1', 
                                        'long_name' : 'Net Phyto Oxygen', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat deptht'},
                            'desc' : 'Net phyto impact on Oxygen (ie. Primary production - repsiration + nitrate reduction)',
                            'f' : lambda x : (x.NPP + x.PhytoNitrateReduction)},
            'ZooResp'   :  { 'req' : ['TotalRespiration_Zoo','TotalRespiration_Gel'] , 
                            'attrs' : {'units'     : 'mmol O2 m-3 s-1', 
                                        'long_name' : 'Zooplankton respiration', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat deptht'},
                            'desc' : 'Total Zoo and Gel respiration',
                            'f' : lambda x : (x.TotalRespiration_Gel + x.TotalRespiration_Zoo)},
            'DOC'    :  { 'req' : ['DCL','DCS'] , 
                            'attrs' : {'units'     : '-', 
                                        'long_name' : 'Dissolved Organic Carbon', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat '},
                            'desc' : 'sum of labile and semi-labie dissolved organic matter components',
                            'f' : lambda x : x.DCL + x.DCS},
            'ZooRespI'    :  { 'req' : ['deptht','ZooResp'] , 
                            'attrs' : {'units'     : '-', 
                                        'long_name' : 'Total Zooplankton respiration - Vertical integral', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat '},
                            'desc' : 'Zooplankton respiration (Gel + Zoo)- Vertical integral',
                            'f' : lambda x : integratevar(x,'ZooResp')},
            'NPPOI'    :  { 'req' : ['deptht','NPPO'] , 
                            'attrs' : {'units'     : '-', 
                                        'long_name' : 'Total Zooplankton respiration - Vertical integral', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat '},
                            'desc' : 'Net phyto impact on Oxygen - Vertical integral',
                            'f' : lambda x : integratevar(x,'NPPO')},
            'OXIDATIONBYDOXI'    :  { 'req' : ['deptht','OXIDATIONBYDOX'] , 
                            'attrs' : {'units'     : '-', 
                                        'long_name' : 'OXIDATIONBYDOX - Vertical integral', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat '},
                            'desc' : 'OXIDATIONBYDOX (nitrifcation and ODU oxidation) - Vertical integral',
                            'f' : lambda x : integratevar(x,'OXIDATIONBYDOX')},
            'bac_oxygenconsumptionI'    :  { 'req' : ['deptht','bac_oxygenconsumption'] , 
                            'attrs' : {'units'     : '-', 
                                        'long_name' : 'Bacterial Oxygen Consumption - Vertical integral', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat '},
                            'desc' : 'Bacterial Oxygen Consumption - Vertical integral',
                            'f' : lambda x : integratevar(x,'bac_oxygenconsumption')},
            'CHLI'    :  { 'req' : ['CHL','deptht'] , 
                            'attrs' : {'units'     : 'mg Chl m-2', 
                                        'long_name' : 'Chlorophyll - vertically integrated', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Chlorophyll - vertically integrated',
                            'f' : lambda x : integratevar(x,'CHL')},
            'VOID'    :  { 'req' : [] , 
                            'attrs' : {'units'     : '', 
                                        'long_name' : '', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat deptht'},
                            'desc' : '',
                            'f' : lambda x : x}
}


'''

etc... Couldn't we fid a way to make this as a recursive function for any type of variable ? 
r2.dom_c.data=r2.dom_c.data*1e6
r2.dom_c.attrs["units"] = 'mmol C m-3'
'''

def add2D(x,keys, verbose=True):
    """
    Add the diagnostic KEY to the input xarray x

    Args:
        x (xarray): xarray containing model outputs
        keys (string or vector of strings): identifier of the diagnostic. There should be a corresponding entry in the diagnsotic dictionnary. 

    Returns:
        xarray : xarray completed with the the diagnostic key
    """    
    if type(keys) is not list:
        keys=[keys]

    for key in keys:
        # Test depedencies
        for d in ddiag2D[key]['req']:
            if d in x.keys():
                continue
            else:
                if verbose:print( 'Lacking ' + d + ' to compute '+ key)
                x= add2D(x,d, verbose=verbose)
        x[key]  = ddiag2D[key]['f'](x) 
        x[key].attrs=ddiag2D[key]['attrs']
        if verbose:print ('just added '+  key +' :' + ddiag2D[key]['desc'])
    return x


def integratevar(x,v,upper=None, lower=None):
    from xgcm import Grid

    grid = Grid(x, 
            coords={"Z": {"center": "deptht", "outer": x['deptht'].attrs['bounds']}},
            metrics = {('Z',):['h']})

    if lower is not None:
        return  grid.integrate(x[v].where(x.deptht<lower),'Z')
    elif upper is not None:
        return  grid.integrate(x[v].where(x.deptht>upper),'Z')
    elif (upper is not None) and (lower is not None):
        return  grid.integrate(x[v].where((x.deptht<lower)&(x.deptht>upper)),'Z')
    else:
        return  grid.integrate(x[v],'Z')


def extentwhere(x,v,condition, value):
    from xgcm import Grid

    grid = Grid(x, 
            coords={"Z": {"center": "deptht", "outer": x['deptht'].attrs['bounds']}},
            metrics = {('Z',):['h']})

    if condition=='lower':
        x['cond']=x[v]<value
        return  grid.integrate(x['cond'],'Z')
    else:
        print('condition unknwon') 

def averagevar(x,v,upper=None, lower=None, conditions=None):
    from xgcm import Grid

    grid = Grid(x, 
            coords={"Z": {"center": "deptht", "outer": x['deptht'].attrs['bounds']}},
            metrics = {('Z',):['h']})    
    if conditions is None:
        X=x[v]
    else:
        X=x[v].where(conditions)

    if lower is not None:
        return  grid.average(x[v].where(x.deptht<lower),'Z')
    elif upper is not None:
        return  grid.average(x[v].where(x.deptht>upper),'Z')
    elif ((upper is not None) & (lower is not None)):
        return  grid.average(X.where((x.deptht<lower) & (x.deptht>upper)),'Z')
    else:
        return  grid.average(X,'Z')

def derivate(x,v,upper=None, lower=None):
    from xgcm import Grid

    grid = Grid(x, 
            coords={"Z": {"center": "deptht", "outer": x['deptht'].attrs['bounds']}},
            metrics = {('Z',):['h']})

    xouter=grid.interp(x[v],'Z')

    return  grid.derivative(xouter,'Z', boundary='extend')

def diaglist(keys=ddiag2D.keys()):
    for k in keys:
        print( "{0:<10}".format(k) + ' - [' + ddiag2D[k]['attrs']['units'] + '] : ' + ddiag2D[k]['desc'])
        print( '  Requires : ' + " ; ".join(ddiag2D[k]['req']))
        # for r in ddiag2D[k]['req']:
        #     print( '          '+r )
        print('\n')
