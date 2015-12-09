from __future__ import print_function, division

import os
import os.path as osp
import matplotlib.pyplot as plt

from warnings import warn
import datetime, time
import glob as gb
from six import string_types
import argparse
import json
import time
import re

import numpy as np
import scipy.linalg as lin
import scipy.stats as sst

import smpce_data_to_corr as smp
import layout as lo
import utils._utils as ucr
import correlation2results as c2r


def _get_conditions_patterns():
    conds_pat = ['NONE.npz', 'LOW.npz', 'MED.npz', 'HIGH.npz']
    return conds_pat

def _get_conditions():
    """
    return ordered conditions
    """
    return ['none', 'low', 'med', 'high']    
    
def ordered_conds():
    return _get_conditions()

def get_signals_filenames(dbase, addjson=None, condpat='', dstate={'sub':1}):
    """
    TODO : loop over subjects ? not sure as our analyses will be subject per subject

    return sorted filenames for addjson (a specific analysis) and condition condpat
    condpat : can be none, NONE, none.npz, or NONE.npz (or any other : low, med, high) 
    """

    params = smp.get_params(dbase, addjson=addjson)

    pth = lo._get_pth(params['layout'], 'signals', glob=True)
    glb = lo._get_glb(params['layout'], 'signals', glob=True)
    
    result_files = gb.glob(osp.join(pth,glb))
    if not result_files:
        raise ValueError("\npth {} \nglb {} \nare not giving anything".format(pth,glb)) 

    if condpat:
        condname, cond_extension = osp.splitext(condpat)
        if not cond_extension: cond_extension = '.npz'
        condpat = condname.upper() + cond_extension
        assert condpat in _get_conditions_patterns()
        return [f for f in sorted(result_files) if f.endswith(condpat)]
    else: 
        return result_files

def get_conds_filenames(dbase, addjson=None, dstate={'sub':1}, pipeline_dic={}):
    """
    TODO : loop over subjects ?
    returns:
    ---------
    conds: dict
        conds['none'] is the list of files of the 'none' condition, etc
    pipeline: string
        taken from the addjson: will take this name, or "default" if addjson=None
    """
    conds = {}
    pipeline_name = 'default'
    if addjson: pipeline_name, ext = osp.splitext(addjson)
    if pipeline_name in pipeline_dic:
        pipeline_name = pipeline_dic[pipeline_name]
    for co in _get_conditions():
        conds[co] = get_signals_filenames(dbase, addjson=addjson, condpat=co, dstate=dstate)
    return conds, pipeline_name




def test_get_npz_filenames(dbase):
    params = smp.get_params(dbase, addjson=None)
    assert params['layout']['signals']['pre'] == 'signal_'
    assert params['analysis']['compute_sess_mask']

    params = smp.get_params(dbase, addjson='NR_0')
    assert params['layout']['signals']['pre'] == 'NR_0_signal_'
    assert not params['analysis']['compute_sess_mask']

    params = smp.get_params(dbase, addjson='NR_1')
    assert params['layout']['signals']['pre'] == 'NR_1_signal_'
    assert not params['analysis']['compute_sess_mask']

    # dstate = {'sub':1}
    pth = lo._get_pth(params['layout'], 'signals', glob=True)
    assert pth == osp.join(dbase, 'sub*/sess*/preproc/extracted_signals')
    glb = lo._get_glb(params['layout'], 'signals', glob=True)
    assert glb == 'NR_1_signal_sub*_sess*_run*'

    return True

def test_get_conds_filenames(dbase, pipelines, nsess=13, ncond=4):
    """
    This function checks for each pipeline:
        1. that the number of files returned per condition is the number of session
        2. that all these files have different content 
    """

    def sha1_of_file(filepath):
        from hashlib import sha1
        with open(filepath, 'rb') as f:
            return sha1(f.read()).hexdigest()
    
    shas = []
    npipe = len(pipelines)
    for pipe in pipelines:
        # print(pipe, len(shas))
        conds, pipeline = get_conds_filenames(dbase, addjson=pipe)
        # print("\n".join(conds['none']))
        check_nb_files = [np.asarray(conds[cond]).shape[0] for cond in conds]
        # we should have nsess files per condition
        assert check_nb_files == [nsess]*ncond
        
        for cond in conds:
            shas.extend([sha1_of_file(fn) for fn in conds[cond]])
            
    # print(shas)
    assert len(set(shas)) == nsess*ncond*npipe, "{} {}".format(len(set(shas)), nsess*ncond*npipe)
    
    return True

def _get_common_labels(conds, idx0=0):
    """
    Returns the common labels across sessions 
    For one session, the conditions should have the same regions, hence
    it is enough to get the regions common across sessions for one condition.
    idx0 provides with the index of that condition. 
    
    Parameters:
    -----------
    conds: dict
        a dictionary with keys the different conditions (eg 'none', 'low', ...)
        each key contains number of sessions npz filenames like :
        "signal_run??_`cond`.npz"
    idx0: int
        index of the condition to determine the common rois across sessions 
    returns:
    --------
    """
    
    cond0 = conds.keys()[idx0] 
    nb_sess = len(conds[cond0])
    
    lsets = []
    for sess in range(nb_sess):
        lsets.append( set((np.load(conds[cond0][sess]))['labels_sig']) )

    return set.intersection(*lsets)

def _test_common_labels(conds):
    """
    check that all conditions have the same common rois (ie labels)
    4 conditions
    """
    common_labels = _get_common_labels(conds, idx0=0)
    assert common_labels == _get_common_labels(conds, idx0=1)
    assert common_labels == _get_common_labels(conds, idx0=2)
    assert common_labels == _get_common_labels(conds, idx0=3)
    return True 

def compute_corr_mtx(conds, common_labels):
    """
    returns
    -------
    arr: list of 4 np.array, each (nsess, common_len, common_len)
    stored_params: dict
        parameters found in the extracted signals 
    """
    conditions_ = conds.keys()
    conditions  = _get_conditions() # ordered
    assert set(conditions_) == set(conditions)
    
    cond0, sess0 = conditions[0], 0 # cond0 should be "none"
    nb_sess = len(conds[cond0])
    shape_c = (nb_sess, len(common_labels), len(common_labels))
    
    conds_arr = {}
    
    for cond in conditions:
        arr_c = np.empty(shape_c)
        arr_c.fill(10.)
        
        for sess in range(nb_sess):
            dctsig = np.load(conds[cond][sess])
            # create a boolean array indicating which of the labels in dctsig are in
            # common_labels
            idx_lab = np.asarray([lab in common_labels 
                                      for lab in dctsig['labels_sig']])
            assert idx_lab.sum() == len(common_labels), "\n {}, {}".format(
                                           idx_lab.shape[0] ,len(common_labels))
            # compute correlation matrix
            arr_c[sess] = np.corrcoef(dctsig['arr_sig_f'][:,idx_lab].T)
            # check all values are ok
            if arr_c[sess].max() > 1.0:
                print('sess {} cond {} '.format(sess,cond))
                print("max should be <= 1.0 {}".format(arr_c[sess].max()))
                print(idx_lab.sum())
                print(np.where(arr_c[sess] > 1.))
                raise ValueError

        # redundant with raise in the loop
        # assert arr_c.max() <= 1.0, "max should be <= 1.0 and is: {}".format(arr_c.max())  
        conds_arr[cond] = arr_c
    
    conds_arr['labels'] = common_labels
    stored_params = np.load(conds[cond0][sess0])['params']
            # should be identical now across cond and sess, 
            # check in the future if this will be overloaded
    
    return conds_arr, stored_params


# functions that work with the dictionary conds
#-------------------------------------------------
# def smpce_mean_cond(cond_arr, cond):
#     """ 
#     Estimates the bias of each condition
#     assumes that cond_arr['cond'] is (nb_of_sess, nb_roi, nb_roi)
#     """
#     # two last dimension must be indentical
#     assert cond_arr[cond].shape[-1] == cond_arr[cond].shape[-2]  
#     return cond_arr[cond].mean(axis=0)

    
def fisher_transf(rho):
    """ take a coefficient of correlation and z transform it 
        see en.wikipedia.org/wiki/Fisher_transformation
    """
    assert rho.max() < 1. and rho.min() > -1. , \
           "rho.max(), rho.min() {} {}".format(rho.max(), rho.min())    
    return (0.5 * np.log((1. + rho) / (1. - rho)))


def summary_per_corr_mtx(cond_arr, summary=np.mean, fisher_transform=True):
    """ 
    input: 
    -----
    array of correlations (nsess, nroi, nroi)
    or (nroi, nroi)
    output:
    -------
    array (nsess,)
    """
    
    if cond_arr.ndim == 2:
        nsess = 1
        nroi = cond_arr.shape[0]
        assert cond_arr.shape[0] == cond_arr.shape[1]
        cond_arr = cond_arr.reshape((1,nroi,nroi))

    elif cond_arr.ndim == 3:
        nsess = cond_arr.shape[0]
        nroi = cond_arr.shape[1]
        assert cond_arr.shape[1] == cond_arr.shape[2]

    iupper = np.triu_indices(nroi,1)

    cval = np.zeros((nsess,)) 
    for sess in range(nsess):
        sess_arr = cond_arr[sess][iupper]
        if fisher_transform:
            sess_arr = fisher_transf(sess_arr)
        cval[sess] = summary(sess_arr)

    return cval


def fisher_transf_cond_arr(cond_arr):
    """
    cond_arr: array
        a (nsess, nroi, nroi) array
    returns:
    --------
        a (nsess, npairs) array    
    """
    assert cond_arr.ndim == 3
    nsess = cond_arr.shape[0]
    nroi = cond_arr.shape[1]
    assert cond_arr.shape[1] == cond_arr.shape[2]

    iupper = np.triu_indices(nroi,1)
    zcond_arr = np.zeros((nsess, len(iupper[0]))) 
    assert len(iupper[0]) == float(nroi*(nroi-1))/2.

    for sess in range(nsess):
        sess_arr = cond_arr[sess][iupper]
        zcond_arr[sess] = fisher_transf(sess_arr)

    return zcond_arr

def conds_2_zconds(conds_arr, ordered_conds, verbose=False):
    """
    """
    zconds_arr = {}
    for ke in ordered_conds:
        if verbose: print("in: ", ke, conds_arr[ke].shape)
        zconds_arr[ke] = fisher_transf_cond_arr(conds_arr[ke])
        if verbose: print("out: ", ke, zconds_arr[ke].shape)
    
    return zconds_arr

def de_bias(zconds_arr, ordered_conds):
    """ 
    Estimates the bias of each condition
    assumes that cond_arr['cond'] is (nb_of_sess, nb_pairs_roi)
    input:
    ------
    zconds: dict
        each key 'cond' contains a (nsess, npairs) numpy array
    returns:
    --------
        zconds array with mean (across sessions) of the 'none' 
        condition removed
    """
    
    assert ordered_conds[0] == 'none'
    estTrueC = zconds_arr['none'].mean(axis=0)
    un_biased = {}
    for ke in ordered_conds:
        un_biased[ke] = zconds_arr[ke] - estTrueC
    
    return un_biased


def remove_cond(zconds_arr, ordered_conds, remove='none'):
    """
    """
    minus_remove = {}
    zremove = zconds_arr[remove] 
    for ke in ordered_conds:
        minus_remove[ke] = zconds_arr[ke] - zremove

    return minus_remove



def summary_per_zconds(zconds_arr, summary=np.mean, verbose=False, ord_conds=ordered_conds()):
    """ 
    summary across npairs dimension (axis=1)
    input: 
    -----
    dict 
        conds: contains (nsess, npairs) arrays
    output:
    -------
    dict 
        conds: array (npairs,)
    """
    
    if verbose: print(type(zconds_arr))
    if verbose: print(zconds_arr.keys())
    summary_zc = {}
    for ke in ord_conds:
        if verbose: print(ke, zconds_arr[ke].shape, 'axis:',1)
        summary_zc[ke] = summary(zconds_arr[ke])

    return summary_zc


def summary_correlation(zconds_arr, ord_conds=ordered_conds()):
    """
    compute the correlation between conditions and none
    """
    summary_corr = {}
    none_cond_values = zconds_arr['none']
    nsess, npairs = none_cond_values.shape

    for ke in ord_conds:
        lcorr = []
        for sess in range(nsess):
            corr = np.corrcoef(none_cond_values[sess], zconds_arr[ke][sess])
            lcorr.append(corr[0,1])
            
        summary_corr[ke] = np.asarray(lcorr)

    return summary_corr


def smpce_bias(conds_arr, ordered_conds):
    """ 
    Estimates the bias of each condition
    assumes that cond_arr['cond'] is (nb_of_sess, nb_roi, nb_roi)
    """
    
    assert ordered_conds[0] == 'none'
    estTrueC = conds_arr['none'].mean(axis=0)
    bias = {}
    for ke in ordered_conds:
        bias[ke] = conds_arr[ke].mean(axis=0) - estTrueC
    
    return bias

def smpce_bias_std(conds_arr, ordered_conds):
    """ 
    Estimates the bias of each condition
    assumes that cond_arr['cond'] is (nb_of_sess, nb_roi, nb_roi)
    """
    
    assert ordered_conds[0] == 'none'
    estTrueC = conds_arr['none'].mean(axis=0)
    bias = {}
    std = {}
    for ke in ordered_conds:
        bias[ke] = conds_arr[ke].mean(axis=0) - estTrueC
        std[ke] = conds_arr[ke].std(axis=0)
    
    return bias, std

def smpce_corr_btw_sess_cond(cond_arr):
    """
    return a (nsess,nsess) correlation matrix
    cond_arr should be (nsess, nroi, nroi) array
    """
    nsess = cond_arr.shape[0]
    nroi = cond_arr.shape[1]
    iupper = np.triu_indices(nroi,1)
    #print(len(iupper[0]))
    upper_corr = np.zeros((nsess, len(iupper[0])))
    for idx,cor in enumerate(cond_arr):
        upper_corr[idx] = cor[iupper]
        #print(cor[iupper].shape)

    return np.corrcoef(upper_corr)


def smpce_var(cond_arr, ordered_conds):
    """ 
    Estimates the variance of each condition
    assumes that cond_arr['cond'] is (nb_of_sess, nb_roi, nb_roi)
    """
    
    assert ordered_conds[0] == 'none'
    
    means = {}
    vari = {}
    for ke in ordered_conds:
        means[ke] = cond_arr[ke].mean(axis=0)
        vari[ke] = np.sqrt(cond_arr[ke].var(axis=0))

#    for ke in ordered_conds:
#        means[ke] = smpce_mean_cond(cond_arr, ke)
#    for ke in ordered_conds:
#        Cm = cond_arr[ke]
#        vari[ke] = (Cm**2 - means[ke]).mean(axis=0)
        
    return vari
    

# plotting functions that work with the dictionary conds
#--------------------------------------------------------

def plot_pipeline_summary(conds_summary, title, pipeline=None, minall=-.4, maxall=.4):

    f, axes = plt.subplots(1, 4, sharey=True, figsize=(16,4))
    f.subplots_adjust(wspace=0.1)
    f.subplots_adjust(right=0.85)
    titlestr = 'Pipeline ' + pipeline + ' - ' + title
    f.suptitle(titlestr, fontsize=24, fontweight='bold', x=.5, y=.01)
    left, bottom, width, height = .87, 0.12, 0.02, 0.77 
    cbar_ax = f.add_axes([left, bottom, width, height])

    minall_ = np.min(np.asarray([conds_summary[k].min() for k in ordered_conds]))
    maxall_ = np.max(np.asarray([conds_summary[k].max() for k in ordered_conds]))
    if not minall: minall=minall_
    if not maxall: maxall=maxall_

    for axe,k in zip(axes,ordered_conds): 
        m = axe.imshow(conds_summary[k], interpolation='nearest', vmin=minall, vmax=maxall)
        axe.set_title(k, fontsize=20)

    cm = f.colorbar(m, cax=cbar_ax)
    return minall_, maxall_
    


#- def save_results(basedir, analysis_label, params, verbose=False):
#-     """
#-     take basedir and analysis label to store the correlation matrix
#-     and the parameters in a directory determine by 
#-     basedir + params['layout']['res']['dir'] + analysis_label+'_'+timestr
#-     
#-     Parameters:
#-     -----------
#-     basedir: string
#-     analysis_label: string
#-     params: dict
#-     """
#- 
#-     return False
#-     permission = 0o770 # "rwxrws---"
#-     
#-     conds =  get_conds_filenames(basedir, params)
#-     common_labels = _get_common_labels(conds)
#-     conds_arr, stored_params = compute_corr_mtx(conds, common_labels)
#- 
#-     params    = smp.get_params(basedir)
#-     resultdir = params['layout']['res']['dir']
#-     matrixdir = params['layout']['res']['corr']
#-     cpsigdir  = params['layout']['res']['sig']
#-     timestr   = time.strftime("%m%d_%H%M%S")
#-     
#-     res_dir   = osp.join(basedir, resultdir)
#-     label_dir = osp.join(res_dir, analysis_label+'_'+timestr)
#-     matrixdir = osp.join(label_dir, matrixdir)
#-     cpsigdir  = osp.join(label_dir, cpsigdir)
#- 
#-     # make directories - will look like results/"analysis_label"_timestr/...
#-     if not osp.isdir(res_dir):
#-         os.makedirs(res_dir, permission)
#-     os.makedirs(label_dir, permission)
#-     os.makedirs(matrixdir, permission)
#-     os.makedirs(cpsigdir, permission)
#- 
#-     # cp files containing signals to cpsigdir
#-     for cond in conds:
#-         for fn in conds[cond]:
#-             shu.copy(fn, cpsigdir)
#-     
#-     # save correlation matrices
#-     fn_save = osp.join(matrixdir, analysis_label)
#-     np.savez(fn_save, conds=conds, common_labels=common_labels, 
#-             conds_arr=conds_arr, stored_params=stored_params)
#- 
#-     return fn_save



#-------------------------------------------------------------------------------
# main
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    #
    parser = argparse.ArgumentParser()

    # Positional required arguments
    help_base_dir = "The directory containing sub01/sess??/... and json files: \n" + \
                    "analysis_parameters.json  data_parameters.json  directory_layout.json"
    parser.add_argument("base_directory", help=help_base_dir)
    help_analysis_label = "The analysis_label for the results "
    parser.add_argument("analysis_label", nargs='?', help=help_analysis_label, default="nofuture")

    # Optional arguments
    parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
    args = parser.parse_args()
    base_dir = args.base_directory
    analysis_label = args.analysis_label
    verbose = args.verbose
    
    assert osp.isdir(base_dir), '{} not a directory'.format(base_dir)
    if verbose: print("launching analyses on :", base_dir)

    params = get_params(base_dir)

    if verbose: 
        nb_sess = params['data']['nb_sess']
        nb_sub = params['data']['nb_sub']
        nb_run = params['data']['nb_run']
        print("nb_sub: {} nb_sess: {}, nb_run: {}".format(nb_sub, nb_sess, nb_run))

    conds = get_conds_filenames(base_dir)
    if verbose: print("conditions in conds: {}".format(conds.keys()))

    # some checks
    assert _test_common_labels(conds)

    fn_saved = save_results(base_dir, analysis_label, params, verbose=verbose)

    if verbose: print(fn_saved)

