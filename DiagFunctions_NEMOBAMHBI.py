from numpy import append
import xarray as xr
import xgcm

# Catalogue of diagnostics for PISCES outputs
# A. Capet - acapet@uliege.be - Feb 2022

ddiag2D = { 'poc'    :  { 'req' : ['POC','GOC'] , 
                            'attrs' : {'units'     : 'mmol C m-3', 
                                        'long_name' : 'Particulate Organic Carbon', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat deptht'},
                            'desc' : 'Total (dead) organic matter carbon content',
                            'f' : lambda x : (x.POC + x.GOC)},
            'TrophicEfficiency'    :  { 'req' : ['ZOO','ZOO2','PHY','PHY2'] , 
                            'attrs' : {'units'     : '-', 
                                        'long_name' : 'Trophic Efficiency', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat deptht'},
                            'desc' : 'Ratio ZOO/PHYTO',
                            'f' : lambda x : (x.ZOO2 + x.ZOO)/(x.PHY + x.PHY2)},
            'TrophicEfficiencyI'    :  { 'req' : ['deptht','ZOO','ZOO2','PHY','PHY2'] , 
                            'attrs' : {'units'     : '-', 
                                        'long_name' : 'Trophic Efficiency - Vert Int.', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat '},
                            'desc' : 'Ratio Int_{z=0-200}(TotZoo)/Int_{z=0-200}(TotPhyto)',
                            # AC 22032022 - replace lower deptht bounds from negative to positive values. 
                            # TODO : implement consistency check with deptht convention
                            'f' : lambda x : (integratevar(x,'ZOO', lower=200)+ integratevar(x,'ZOO2', lower=200))/(integratevar(x,'PHY', lower=200)+ integratevar(x,'PHY2', lower=200))},
            'ratioLargeM'    :  { 'req' : ['ratioLarge'] , 
                            'attrs' : {'units'     : '-', 
                                        'long_name' : 'Large ratio - Average of the ratio where TotPHY > 0.01 µM ', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Average ( DIA/PHYTO ) where TotPhyto>0.01',
                            'f' : lambda x : averagevar(x,'ratioLarge', conditions=(x.PHY2+x.PHY)>0.01)},                            
            'pocM200'    :  { 'req' : ['poc'] , 
                            'attrs' : {'units'     : 'mmol C m-3', 
                                        'long_name' : 'Average poc content, 0-200m', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Average poc content, 0-200m',
                            'f' : lambda x : averagevar(x,'poc', lower=200)},                            
            'CHLM5'    :  { 'req' : ['CHL'] , 
                            'attrs' : {'units'     : 'mg Chl m-3', 
                                        'long_name' : 'Average Chl content, 0-5m', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Average chl content, 0-5m',
                            'f' : lambda x : averagevar(x,'CHL', lower=5)},                            
            'ratioLarge'    :  { 'req' : ['PHY','PHY2'] , 
                            'attrs' : {'units'     : '-', 
                                        'long_name' : 'Large ratio - DIA/PHYTO ', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat deptht'},
                            'desc' : 'DIA/PHYTO',
                            'f' : lambda x : x.PHY2/(x.PHY2+x.PHY)},                            
            'pocF'    :  { 'req' : ['POC','GOC'] , 
                            'attrs' : {'units'     : 'mg C m-2 d-1', 
                                        'long_name' : 'Downward Flux of Particulate Organic Carbon', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat deptt'},
                            'desc' : 'Downward Flux of Particulate Organic Carbon ! ASSUMED VALUES FOR SINKING VELOCITIES : 2 and 50 m/d !',
                            'f' : lambda x : 12*(x.POC*2 + x.GOC*50)},  # TODO : enable parameter read in fabm.yaml
            'pocF100'    :  { 'req' : ['pocF','deptht'] , 
                            'attrs' : {'units'     : 'mg C m-2 d-1', 
                                        'long_name' : 'Downward Flux of Particulate Organic Carbon @ 100m', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Downward Flux of Particulate Organic Carbon interpolated @ 100m',
                            'f' : lambda x : x.pocF.interp(deptht=100)}, 
            'zchlmax'    :  { 'req' : ['CHL','deptht'] , 
                            'attrs' : {'units'     : 'm', 
                                        'long_name' : 'Depth of Max. Chl.', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Depth of Max Chl (straight computation with idxmax)',
                            'f' : lambda x : x.CHL.idxmax('deptht')},
            'CHLatMAX'    :  { 'req' : ['CHL','deptht'] , 
                            'attrs' : {'units'     : 'mg Chl m-3', 
                                        'long_name' : 'Max. Chl.', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Max Chl along the vertical',
                            'f' : lambda x : x.CHL.max('deptht')},
            'nitracline'    :  { 'req' : ['NO3','deptht'] , 
                            'attrs' : {'units'     : 'm', 
                                        'long_name' : 'Nitracline Depth', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Depth of max nitrcaline gradient',
                            'f' : lambda x : derivate(x,'NO3').idxmax('deptht')},
            # 'MLD'    :  { 'req' : ['rho'] , 
            #                 'attrs' : {'units'     : 'm', 
            #                             'long_name' : 'Mixed Layer Depth', 
            #                             'valid_min' : -1e20,
            #                             'valid_max' : 1e20,
            #                             'cell_methods' : 'time: mean',
            #                             'coordinates': 'lon lat'},
            #                 'desc' : ' First Depth where RHO > RHO(-3m) - 0.125',
            #                 'f' : lambda x : (abs(x.rho - x.rho.interp(deptht=-3)-0.125)).idxmin('deptht')},
            'TPPI'    :  { 'req' : ['TPP','deptht'] , 
                            'attrs' : {'units'     : 'mg C m-2 d-1', 
                                        'long_name' : 'Net Primary Production - vertically integrated', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Net Primary Production - vertically integrated',
                            'f' : lambda x : integratevar(x,'TPP')},
            'CHLI'    :  { 'req' : ['CHL','deptht'] , 
                            'attrs' : {'units'     : 'mg Chl m-2', 
                                        'long_name' : 'Chlorophyll - vertically integrated', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Chlorophyll - vertically integrated',
                            'f' : lambda x : integratevar(x,'CHL')},
            'OMZextent'    :  { 'req' : ['O2','deptht'] , 
                            'attrs' : {'units'     : 'm', 
                                        'long_name' : 'OMZ extent', 
                                        'valid_min' : -1e20,
                                        'valid_max' : 1e20,
                                        'cell_methods' : 'time: mean',
                                        'coordinates': 'lon lat'},
                            'desc' : 'Vertical extent where oxygen <90µM',
                            'f' : lambda x : extentwhere(x,'O2','lower',120)},            # 'CHLsurf'    :  { 'req' : ['total_chlorophyll_calculator_result_'] , 
            #                 'attrs' : {'units'     : 'mg C m-3', 
            #                             'long_name' : 'Surface Chlorophyll', 
            #                             'valid_min' : -1e20,
            #                             'valid_max' : 1e20,
            #                             'cell_methods' : 'time: mean',
            #                             'coordinates': 'lon lat'},
            #                 'desc' : '',
            #                 'f' : lambda x : x.total_chlorophyll_calculator_result.interp(deptht=-1)},
            # 'TotNI_'    :  { 'req' : ['total_nitrogen_calculator_result'] , 
            #                 'attrs' : {'units'     : 'mmol m-2', 
            #                             'long_name' : 'Tot N - vertically integrated', 
            #                             'valid_min' : -1e20,
            #                             'valid_max' : 1e20,
            #                             'cell_methods' : 'time: mean',
            #                             'coordinates': 'lon lat'},
            #                 'desc' : '',
            #                 'f' : lambda x : integratevar(x,'total_nitrogen_calculator_result')},
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
