from __future__ import print_function, division
import os 
import os.path as osp
import json
import re
import glob as gb
from six import string_types

def _get_key_dict(dlayo, entities = ["subjects", "sessions", "runs"]):
    """
    ["subjects", "sessions", "runs"] are the entities of the directory layout
    return a dict with the keys 
    eg: 
    dkeys["subjects"] == 'sub'
    """

    dkeys = {}
    for ent in entities:
        assert ent in dlayo
        assert "key" in dlayo[ent]
        dkeys[ent] = dlayo[ent]["key"]
    
    return dkeys 

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
        if 'base' in dlayo:
            dlayo['base'] = dbase

    return dlayo 

# !!! to be implemented
# sub_str = lo._format(sub_pth, dstate)

def _format(toformat, fdict, failsifnotcomplete=True):
    """
    takes a string to format, a dict that may contain too many keys,
    remove the uncessary keys and returns the formated string. 
    """
    #assert type(toformat) in [str, unicode]
    assert isinstance(toformat, string_types)

    if failsifnotcomplete:
        # check that all keys in the string have a corresponding value in fdic
        #!!! to be implemented
        pass
    return toformat.format(**fdict)

#    reg = '\{([a-z:]*)\}' #should capture the {aa:bb} in the string
#    p = re.compile(reg)
#    args = p.findall(toformat)
#    reduced_fdic = {}
#    for arg in args:
#        if arg:
#            k,v = arg.split(':')
#            assert k in fdict
#            reduced_fdic[k] = fdict[k]
#    

def _get_pth(dlayo, key, verbose=False):
    """
    construct the path, the file name, and the dictionary of keys for the path  
    returns
    -------
    strpth: string
    pthdic: dict 
    """

    try:
        assert key in dlayo
    except:
        print("Key {} not in dict keys {}".format(key, dlayo.keys()))

    strpth = dlayo['base'] 
    pdic = {}
    # first construct the *previous* path of this key with "base" and "key-val"
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
    construct the path, the file name, and the dictionary of keys for the file name 
    returns:
    --------
    the string corresponding to the glob pattern 
    """
    try:
        assert key in dlayo
    except:
        print("Key {} not in dict keys {}".format(key, dlayo.keys()))

    if 'glb' not in dlayo[key]:
        raise ValueError(" 'glb' is not a key for '{}'".format(key))

    # return string if simply a string
    if type(dlayo[key]['glb']) in (str, unicode):
        return dlayo[key]['glb'], {}

    # construct file name using the 'glb' section of dlayo[key] 
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

def _get_pthglb(dlayo, key, verbose=False):
    """
    """
    pth, pth_dict = _get_pth(dlayo, key, verbose=verbose)
    glb, glb_dict = _get_glb(dlayo, key, verbose=verbose)

    pthglb_dict = pth_dict.copy()
    pthglb_dict.update(glb_dict)

    return osp.join(pth, glb), pthglb_dict

def _get_apthglb(dlayo, key, dstate, verbose=False):
    """
    """
    pthglb, _ = _get_pthglb(dlayo, key, verbose=verbose)
    return osp.join(pthglb).format(**dstate)

def _get_apth(dlayo, key, dstate, verbose=False):
    """
    get a path, format it with dstate
    """
    pth, _ = _get_pth(dlayo, key, verbose=verbose)
    return osp.join(pth).format(**dstate)


def _get_oneof(dlayo, key, verbose=False):
    """
    pattern for one directory 
    """
    pth, pth_dict = _get_pth(dlayo, key, verbose=verbose)
    if 'hasdirs' in dlayo[key] and (dlayo[key]['hasdirs'] is True):
        pth = osp.join(pth, dlayo[key]['key'] + dlayo[key]['val'])
        pth_dict[dlayo[key]['key']] = None
    else:
        # only implemented when there are multiple directories
        raise ValueError("{} has not multiple dirs".format(key))
        
    return pth, pth_dict 

def _get_aoneof(dlayo, key, dstate, verbose=False):
    """
    one directory with dstate to instanciate the sub and sess
    """
    pth, _pth_dict = _get_oneof(dlayo, key, verbose=verbose)

    return pth.format(**dstate)

def _get_alistof(dlayo, key, dstate={}, return_idxs=False, verbose=False):
    """
    for list of directories (ie subjects, sessions)
    """
    apth = _get_apth(dlayo, key, dstate,  verbose=verbose)
    apthglb = _get_apthglb(dlayo, key, dstate, verbose=verbose)

    if 'hasdirs' in dlayo[key] and (dlayo[key]['hasdirs'] is True):
        # !!! TODO : here we should not check hasdirs, but if there are several values
        apth_kv = osp.join(apth, dlayo[key]['key'] + dlayo[key]['val'])
        lglb = len(gb.glob(apthglb))
        vals = range(1,lglb+1)  # TODO : make it more general, this is only 
                                # if the values are integers
        if return_idxs:
            return vals, [apth_kv.format(**{dlayo[key]['key']:val}) for val in vals]
        else:
            return [apth_kv.format(**{dlayo[key]['key']:val}) for val in vals]

    else: 
        assert return_idxs == False # cannot return idxs if not a list of dirs
        return gb.glob(apthglb)
        # raise ValueError("{} has not multiple dirs".format(key))

def _get_aunique(dlayo, key, dstate={}, check_exists=True):

    apthglb = _get_apthglb(dlayo, key, dstate, verbose=False)
    files = gb.glob(apthglb)
    if len(files) != 1:
        raise ValueError(" key '{}' is suppose to yield only one file".format(key))
    else:
        if check_exists:
            if not osp.isfile(files[0]):
                raise ValueError(" file/dir '{}' does not exist".format(files[0]))
    return files[0]


    

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




