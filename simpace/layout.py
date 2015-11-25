from __future__ import print_function, division
import os 
import errno
import os.path as osp
import json
import re
import glob as gb
from six import string_types


# Documentation:

"""
The json file contains a set of keys that are the entities that we wish to retrieve.
List of Json file top level keys, also called entities: 
                        ["subjects", 
                         "sessions", 
                         "runs", 
                         "smoothed", 
                         "mvt6params", 
                         ...]

Each entity will contain a set of (ek,v) pairs. ek can be  
    'base'      : if we want a path, this is the last added bit
    'pth'       : made with a set of entities, get the ek,val,base  
    'hasdirs'   : determine if there is another level to the path
    'key'       : used to construct both path and filenames
    'val'       : used to construct both path and filenames
    'glb'       : list of stuff to construct filenames to glob
                    contains 'pre' 'post' 'lnk' and list of ek
    'pre' 'post' 'lnk' : use to construc glb


A entity will contain one of several files or directory, or values.
There are two main usages:
    1. We want to construct the path that leads to a specific file. To do this, 
    we use the k-v: base, pth, hasdir, val, key. 
    The path is constructed with the entities in pth, for each we add : 
        /   'key' + val + base / 'key' + val + base / ...
    If `glob` is True, then val are replaced by '*' and the result can be use to 
    glob the pattern

    2. We want to construct a pattern to glob some *files*
    This takes the glb list and constructs a series of
        dlayo[k]['key'] + dlayo[k]['val'] + dlayo[k]['lnk']


  "smoothed": {
    "base": "smooth",
    "pth": ["subjects", "sessions", "runs"],
    "hasdirs": false,
    "pre": "sra",
    "post": "????.nii*",
    "glb": ["pre", "subjects", "sessions", "runs", "post"]
  },


"""


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

def _get_pth(dlayo, key, glob=False, verbose=False):
    """
    construct the path to this key and the dictionary of keys for the path  
    if glob = True: returns a path that should be glob
    if glob = False: returns a path and a dict, 
              Use the dict dstate to instanciate

    returns
    -------
    strpth: string
    pthdic: dict 
    """

    try:
        assert key in dlayo
    except:
        err_msg = "Key {} not in dict keys {}".format(key, dlayo.keys())
        raise ValueError, err_msg

    strpth = dlayo['base'] 
    pdic = {}
    # first construct the *previous* path of this key with "base" and "key-val"
    for p in dlayo[key]['pth']:
        base = dlayo[p]['base']
        add_dir = ''
        if 'hasdirs' in dlayo[p] and (dlayo[p]['hasdirs'] is True): 
            if glob:
                val = '*'
            else:
                val = dlayo[p]['val']
            # add the k-v liason in the future 
            add_dir = dlayo[p]['key'] + val
            pdic[dlayo[p]['key']] = None
        if add_dir == '':
            strpth = osp.join(strpth, base)
        else:
            strpth = osp.join(strpth, base, add_dir)
        
    # then add the base of this specific key:
    if dlayo[key]['base'] != '':
        strpth = osp.join(strpth, dlayo[key]['base'])

    if glob: 
        return strpth
    else:
        return strpth, pdic 

def _get_glb(dlayo, key, glob=False, verbose=False):
    """
    This constructs a string from the "glb" key that can be used to glob
    It constructs only the file name(s), but not the path to this / these 
    file(s). It also returns a dictionary of keys to be included
    for the file name 
    parameters:
    -----------
    dlayo: dict
        the dictionary that contains the layout from a json file
    key: string
        the entity we are interested in 
    glob: boolean
        return a string that I can glob ?

    returns:
    --------
    the string to be use to glob the pattern 
    """
    try:
        assert key in dlayo
    except:
        print("Key {} not in dict keys {}".format(key, dlayo.keys()))

    if 'glb' not in dlayo[key]:
        raise ValueError(" 'glb' is not a key for '{}', ".format(key) + 
                         " I cannot _get_glb for that key ")

    # return string if simply a string
    if type(dlayo[key]['glb']) in (str, unicode):
         strglb, fdic = dlayo[key]['glb'], {}
    else:
        # construct file name using the 'glb' section of dlayo[key] 
        strglb = ''
        fdic = {}
        for k in dlayo[key]['glb']:
            if k in dlayo[key]:
                strglb += dlayo[key][k]
            elif k in dlayo:
                if glob:
                    val = '*'
                else:
                    val = dlayo[k]['val']
                # for the last k, we may not want to put the lnk - it is now here
                strglb += dlayo[k]['key'] + val + dlayo[k]['lnk']
                fdic[dlayo[k]['key']] = None
            else:
                raise ValueError, "entity: k {} key {} unknown".format(k, key)

    if glob: 
        return strglb
    else:
        return strglb, fdic 

def _get_pthglb(dlayo, key, glob=False, verbose=False):
    """
    parameters:
    -----------
    dlayo: dict
        the dictionary that contains the layout from a json file
    key: string
        the entity we are interested in 
    glob: boolean
        return a string that I can glob ?

    returns:
    --------


    """
    if glob:
        pth = _get_pth(dlayo, key, glob=glob, verbose=verbose)
        glb = _get_glb(dlayo, key, glob=glob, verbose=verbose)
        return osp.join(pth, glb), None
    else:
        pth, pth_dict = _get_pth(dlayo, key, glob=glob, verbose=verbose)
        glb, glb_dict = _get_glb(dlayo, key, glob=glob, verbose=verbose)
        pthglb_dict = pth_dict.copy()
        pthglb_dict.update(glb_dict)
        return osp.join(pth, glb), pthglb_dict


def _get_apthglb(dlayo, key, dstate, glob=False, verbose=False):
    """
    get the pth and glob, then instanciate this with the dstate

    """
    pthglb, _ = _get_pthglb(dlayo, key, glob=glob, verbose=verbose)
    try:
        apthglb = osp.join(pthglb).format(**dstate)
    except:
       raise ValueError(  " pthglb is {} ".format(pthglb) 
                        + " dstate is {} ".format(dstate)
                        + " and pthglb cannot be formated with dstate.")
    return apthglb

def _get_apth(dlayo, key, dstate, verbose=False):
    """
    get a path, format it with dstate
    """
    pth, _ = _get_pth(dlayo, key, glob=False, verbose=verbose)
    return osp.join(pth).format(**dstate)

def _get_oneof(dlayo, key, verbose=False):
    """
    pattern for one directory 
    first get the path for the entity, then check that this entity has
    directories, last append entity_key, entity_val to this path
    example : /home/../sub{sub:02d} for entity 'subjects'

    returns:
    --------
    pth: string
        the path for one of the directory
    pth_dict: dict
        contains keys to be filled to get an actual path
    """
    pth, pth_dict = _get_pth(dlayo, key, glob=False, verbose=verbose)
    if 'hasdirs' in dlayo[key] and (dlayo[key]['hasdirs'] is True):
        pth = osp.join(pth, dlayo[key]['key'] + dlayo[key]['val'])
        pth_dict[dlayo[key]['key']] = None
    else:
        # only implemented when there are multiple directories
        raise ValueError("{} has not multiple dirs".format(key))
        
    return pth, pth_dict 

def _get_aoneof(dlayo, key, dstate, verbose=False):
    """
    Get oneof and instanciate with dstate 
    """
    pth, _pth_dict = _get_oneof(dlayo, key, verbose=verbose)

    return pth.format(**dstate)

def _get_alistof(dlayo, key, dstate={}, return_idxs=False, verbose=False):
    """
    for list of directories (ie subjects, sessions)
    """
    apth = _get_apth(dlayo, key, dstate,  verbose=verbose)
    apthglb = _get_apthglb(dlayo, key, dstate, glob=False, verbose=verbose)

    if 'hasdirs' in dlayo[key] and (dlayo[key]['hasdirs'] is True):
        # !!! TODO : here we should not check hasdirs, but if there are several values
        apth_kv = osp.join(apth, dlayo[key]['key'] + dlayo[key]['val'])
        lglb = len(gb.glob(apthglb))
        vals = range(1, lglb+1)  # TODO : make it more general, this is only 
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
    """
    first _get_aphtglb then check if only one element (file) is returned
    """
    apthglb = _get_apthglb(dlayo, key, dstate, glob=False, verbose=False)
    files = gb.glob(apthglb)
    if len(files) != 1:
        raise ValueError(" key '{}' is suppose to yield exactly one file".format(key) +
                         "\n and is : {}".format(files) + 
                         "\n with apthglb on pattern {}".format(apthglb))
    else:
        if check_exists:
            if not osp.isfile(files[0]):
                raise ValueError(" file/dir: '{}' does not exist".format(files[0]) +
                                 "\n with glob on pattern {}".format(apthglb))
    return files[0]



#--------------------------------------- a few utilities functions ---------------# 

def _check_dict_full(adic, verbose=False):
    """
    """
    for k in adic:
        if adic[k] is None: return False 
    return True


def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z


def symlink_force(target, link_name):
    try:
        os.symlink(target, link_name)
    except OSError, e:
        if e.errno == errno.EEXIST:
            os.remove(link_name)
            os.symlink(target, link_name)


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



#    
#    
#    
#    

