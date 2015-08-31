from __future__ import print_function, division
import os 
import os.path as osp
import json
import glob as gb

def get_layout(jsonfile, dbase=None, verbose=False):
    """
    parameters:
    -----------
    dbase:  string
            base directory containing subjects directory
    jsonfile:  string
    returns:
    --------
    dict
    """

    # Read json files at the base directory
    if dbase is None:
        dbase = "."

    if not osp.isfile(jsonfile):
        msg = " {} is not a file".format(jsonfile)
        raise ValueError(msg)

    with open(jsonfile) as fjson:
        dlayo = json.load(fjson)

    if not osp.isdir(dbase):
        msg = " {} is not a directory".format(dbase)
        raise ValueError(msg)
    else:
        dlayo['base'] = dbase

    return dlayo 

def _get_pth(dlayo, key, verbose=False):
    """
    construct the path, the file name, and the dictionary of keys for the path and file name 
    """

    try:
        assert key in dlayo
    except:
        print("Key {} not in dict keys {}".format(key, dlayo.keys()))

    strpth = dlayo['base'] 
    pdic = {}
    # first construct the previous path of this key with "base" and "dirs"
    for p in dlayo[key]['pth']:
        base = dlayo[p]['base']
        add_dir = ''
        if 'hasdirs' in dlayo[p] and (dlayo[p]['hasdirs'] is True): 
            # add the k-v liason in the future 
            add_dir = dlayo[p]['key']+dlayo[p]['val']
            pdic[dlayo[p]['key']] = None
        strpth = osp.join(strpth, base, add_dir)
        
    # then add the base of this specific key:
    strpth = osp.join(strpth, dlayo[key]['base'])

    return strpth, pdic 


def _get_glb(dlayo, key, verbose=False):
    """
    construct the path, the file name, and the dictionary of keys for the path and file name 
    """
    try:
        assert key in dlayo
    except:
        print("Key {} not in dict keys {}".format(key, dlayo.keys()))

    # return string if simply a string
    if type(dlayo[key]['glb']) in (str, unicode):
        return dlayo[key]['glb'], {}

    # construct file name using 'glb' 
    strglb = ''
    fdic = {}
    for k in dlayo[key]['glb']:
        if k in dlayo[key]:
            strglb += dlayo[key][k]
        elif k in dlayo:
            strglb += dlayo[k]['key'] + dlayo[k]['val'] + dlayo[k]['lnk']
            fdic[dlayo[k]['key']] = None
        else:
            raise ValueError, "k {} key {} unknown".format(k, key)

    return strglb, fdic

def _check_dict_full(adic, verbose=False):
    """
    """
    for k in adic:
        if adic[k] is None: return False 
    return True





# def get_section(dlayo, section, verbose=False):
#     """
#     returns a section of the json file. Does not work for attribute
#     parameters:
#     -----------
#     dlayo: dict 
#            dict with the directory structure    
#     returns:
#     --------
#     list of string
#     """
#     if verbose:
#         print('>--------------------------------------')
#         print(dlayo)
#     if type(dlayo) is not dict:      # eg ['str','unicode','list']
#         if verbose: print('is not dict')
#         return False
#     if len(dlayo.keys()) == 0:
#         if verbose: print('len(dlayo.keys()) == 0')
#         return False
#     if section in dlayo.keys():
#         if verbose: print('<--- section found')
#         return dlayo[section]
#     else: 
#         for k in dlayo.keys():
#             found = get_section(dlayo[k], section, verbose=verbose)
#             if found is False:
#                 continue
#             else:
#                 if verbose: print('found it',found)
#                 return found
#         if verbose:
#             print('{} not in keys'.format(section))
#             print('<---')
#         return False

# def get_attribute(dlayo, section, attribute, verbose=False):
#     """
#     returns an attribute section of the json file. Does not work for attribute
#     parameters:
#     -----------
#     dlayo: dict 
#            dict with the directory structure    
#     returns:
#     --------
#     list of string
#     """
#     
#     the_section = get_section(dlayo, section, verbose=False)
#     try:
#         return the_section[attribute]
#     except:
#         return None




