
# coding: utf-8

# In[1]:


import marvin
from marvin.tools.maps import Maps

# Returns nested dictionary of the map, inverse variances, and mask of each relevent emission line
def getDataMap(plate,ifu):# plate number and ifu number
    plateifu = plate + '-' + ifu
    cube = marvin.tools.Cube(plateifu)
    
    maps = Maps(plateifu)
    HaF_map = maps["emline_gflux_ha_6564"]
    HbF_map = maps["emline_gflux_hb_4862"]
    OII_map = maps["emline_gflux_oii_3727"]
    OIII_map = maps["emline_gflux_oiii_5008"]
    NII_map = maps["emline_gflux_nii_6585"]
    
    HaF = HaF_map,HaF_map.ivar,HaF_map.mask
    HbF = HbF_map,HbF_map.ivar,HbF_map.mask
    OII = OII_map,OII_map.ivar,OII_map.mask
    OIII = OIII_map,OIII_map.ivar,OIII_map.mask
    NII = NII_map,NII_map.ivar,NII_map.mask
    
    HaF = dict(zip(['map','ivar','mask'],HaF))
    HbF = dict(zip(['map','ivar','mask'],HbF))
    OII = dict(zip(['map','ivar','mask'],OII))
    OIII = dict(zip(['map','ivar','mask'],OIII))
    NII = dict(zip(['map','ivar','mask'],NII))
    
    #     return HaF,HbF,OII,OIII,NII
    
    return{
        'HaF': HaF,
        'HbF': HbF,
        'OII': OII,
        'OIII': OIII,
        'NII': NII
    }
    
# Example of accessing data
getDataMap('8485','1901')['HaF']['mask'][17]

