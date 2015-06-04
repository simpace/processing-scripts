from __future__ import print_function, division

import os
import os.path as osp
import glob as gb
import shutil as shu
import argparse
import time

import numpy as np

from smpce_data_to_corr import get_params


def _get_signals_filenames(basedir, params):
    
    layo = params['layout']
    nb_sess = params['data']['nb_sess']
    #nb_sub = params['data']['nb_sub']
    #nb_run = params['data']['nb_run']
    druns = layo['dir']['runs']
    dsig = layo['dir']['signals']
    
    conditions = ['none', 'low', 'med', 'high']
    conds = {}
    conds_pat = ['*NONE.npz', '*LOW.npz', '*MED.npz', '*HIGH.npz']
    for c in conditions:
        conds[c] = []

    sessions = range(1,nb_sess+1)

    for sess in sessions:
        sesstr = 'sess{:02d}'.format(sess)
        fulldsig = osp.join(basedir, 'sub01', sesstr, druns, dsig)
        for idx, c in enumerate(conditions):
            conds[c].append(gb.glob(osp.join(fulldsig, conds_pat[idx]))[0])
        
    # some basic checks: all the same length
    assert not np.any( np.diff(np.asarray([len(conds[c]) for c in conditions])) )
    
    return conds

def ordered_conds():
    return ['none', 'low', 'med', 'high']

def _get_common_labels(conds, idx0=0):
    
    cond0 = conds.keys()[idx0] 
    nb_sess = len(conds[cond0])
    
    lsets = []
    for sess in range(nb_sess):
        lsets.append( set((np.load(conds[cond0][sess]))['labels_sig']) )

    return set.intersection(*lsets)

def compute_corr_mtx(conds, common_labels):
    """
    returns
    -------
    arr: list of 4 np.array, each (nsess, common_len, common_len)
    """
    idx0 = 0
    conditions_ = conds.keys()
    conditions  = ordered_conds()
    assert set(conditions_) == set(conditions)
    
    cond0 = conditions[idx0]
    sess0 = 0
    nb_sess = len(conds[cond0])
    #nb_cond = len(conds)
    shape_c = (nb_sess, len(common_labels), len(common_labels))
    
    conds_arr = {}
    
    for cond in conditions:
        arr_c = np.empty(shape_c)
        arr_c.fill(10.)
        
        for sess in range(nb_sess):
            dctsig = np.load(conds[cond][sess])
            idx_lab = np.asarray([lab in common_labels for lab in dctsig['labels_sig']])
            # com_lab = [lab for lab in dctsig['labels_sig'] if lab in common_labels]
            arr_c[sess] = np.corrcoef(dctsig['arr_sig_f'][:,idx_lab].T)
            
            assert idx_lab.sum() == len(common_labels), "{},{}".format(
                                idx_lab.shape[0] ,len(common_labels))
            

        conds_arr[cond] = arr_c
    
    conds_arr['labels'] = common_labels
    stored_params = np.load(conds[cond0][sess0])['params']
            # should be identical now across cond and sess check in the future if overloaded
    
    return conds_arr, stored_params


def save_results(basedir, analysis_label, params, verbose=False):

    permission = 0o770
    def ordered_conds():
        return ['none', 'low', 'med', 'high']
    
    conds = _get_signals_filenames(basedir, params)
    common_labels = _get_common_labels(conds)
    conds_arr, stored_params = compute_corr_mtx(conds, common_labels)

    params = get_params(basedir)
    resultdir = params['layout']['res']['dir']
    matrixdir = params['layout']['res']['corr']
    mvsigdir = params['layout']['res']['sig']
    timestr = time.strftime("%m%d_%H%M%S")
    
    res_dir = osp.join(basedir, resultdir)
    if not osp.isdir(res_dir):
        os.makedirs(res_dir, permission)
    
    # make directories - will look like results/"analysis_label"_timestr/...
    label_dir = osp.join(res_dir, analysis_label+'_'+timestr)
    matrixdir = osp.join(label_dir, matrixdir)
    mvsigdir = osp.join(label_dir, mvsigdir)

    os.makedirs(label_dir, permission)
    os.makedirs(matrixdir, permission)
    os.makedirs(mvsigdir, permission)

    # cp files to mvsigdir
    for cond in conds:
        for fn in conds[cond]:
            shu.copy(fn, mvsigdir)
    
    # save correlation matrices
    fn_save = osp.join(matrixdir, analysis_label)
    np.savez(fn_save, conds=conds, common_labels=common_labels, 
            conds_arr=conds_arr, stored_params=stored_params)

    return fn_save

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

    conds = _get_signals_filenames(base_dir, params)

    if verbose: print(conds.keys())
    common_labels = _get_common_labels(conds)
    assert common_labels == _get_common_labels(conds, idx0=3)
    assert common_labels == _get_common_labels(conds, idx0=2)

    fn_saved = save_results(base_dir, analysis_label, params, verbose=verbose)

    if verbose: print(fn_saved)


