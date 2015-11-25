import os.path as osp
import json
import numpy as np
import argparse
#from six import string_types
from nilearn import masking as msk
from nilearn._utils import concat_niimgs
# from nilearn.image.image import _compute_mean
import utils._utils as ucr
import utils.setup_filenames as suf 
import layout as lo

DIRLAYOUT = 'directory_layout.json'
DATAPARAM = 'data_parameters.json'
ANALPARAM = 'analysis_parameters.json'

def get_params(dbase, verbose=False):
    """
    Read json files at the base directory dbase
    
    parameters:
    -----------
    dbase:  string
            base directory containing subjects directory
            *and* jason files
    """

    # Read json files at the base directory
    fn_layo = osp.join(dbase, DIRLAYOUT)
    fn_data = osp.join(dbase, DATAPARAM)
    fn_anal = osp.join(dbase, ANALPARAM)

    dlayo = lo.get_layout(fn_layo, dbase=dbase) #json.load(flayo)
    with open(fn_data) as fdata:
        ddata = json.load(fdata)
    with open(fn_anal) as fanal:
        danal = json.load(fanal)

    params = {'layout': dlayo, 'data': ddata, 'analysis': danal}
    return params


def process_all(dbase, params=None, verbose=False):
    """
    parameters:
    -----------
    dbase:  string
            base directory containing subjects directory
            and jason files
    """
    if not params:
        params = get_params(dbase, verbose=verbose)
    
    dlayo = params['layout']
    # ddata = params['data']

    # make state dict:
    dkeys = lo._get_key_dict(dlayo, entities = ["subjects", "sessions", "runs"])
    dstate = {}
    for k in dkeys: 
        dstate[dkeys[k]] = None
   
    idxs, sub_dirs = lo._get_alistof(dlayo, "subjects", dstate, return_idxs=True)

    subjs_info = {}
    for sub_idx, sub_dirs in zip(idxs, sub_dirs): # start idx at 1
        dstate[dkeys["subjects"]] = sub_idx
        subject_str = lo._get_aoneof(dlayo, "subjects", dstate)
        if verbose: print("subject_str, dstate:",  subject_str, dstate)
        subjs_info[subject_str] = do_one_subject(dstate, dkeys, params, verbose=verbose) 

    return subjs_info

def do_one_subject(dstate, dkeys, params, verbose=False):
    """
    launch sessions processing for current subject 

    parameters:
    -----------
    dstate: dict 
            contains subject index 
    dkeys: dict 
            contains keys for entities  ["subjects", "sessions", "runs"]
    params: dict
            parameters for layout, data and analysis
            
    """
    dlayo = params['layout']
    idxs, sess_dirs = lo._get_alistof(dlayo, "sessions", dstate, return_idxs=True)

    sesss_info, sess_curr = {}, {}
    for sess_idx, sess_dir in zip(idxs, sess_dirs): 
        dstate[dkeys["sessions"]] = sess_idx
        sess_curr['sess_dir'] = sess_dir
        if verbose: print('\n' + '-'*33  + "\n" + sess_dir)
        if verbose: print("sess_dir, dstate:",  sess_dir, dstate)
        sesss_info[sess_dir] = do_one_sess(dstate, dkeys, params, verbose=verbose) 

    return sesss_info
    
def do_one_sess(dstate, dkeys, params, verbose=False):
    """
    launch runs processing for sess_curr 

    parameters:
    -----------
    dstate: dict 
            contains subject and sess index 
    dkeys: dict 
            contains keys for entities  ["subjects", "sessions", "runs"]
    params: dict
            parameters for layout, data and analysis
    """            

    dlayo = params['layout']
    nb_runs = params['data']['nb_run'] 
    assert nb_runs == 4 # 4debug

    # get things that should exist - improve naming scheme ...
    ptr = {}
    ptr['runs_dir'] = lo._get_apth(dlayo, "runs", dstate) # will be preproc
    ptr['smoothed_dir'] = lo._get_apth(dlayo, "smoothed", dstate) # ../smooth
    ptr['aal_dir'] = lo._get_apth(dlayo, "aal_roi", dstate) #  preproc/..
    ptr['aal_files'] = lo._get_alistof(dlayo, "aal_roi", dstate) #  
    ptr['signals_dir'] = lo._get_apth(dlayo, "signals", dstate)
    ptr['mvtfile'] = lo._get_aunique(dlayo, "mvt6params", dstate)
    ptr['csf_file'] = lo._get_aunique(dlayo, "csf_mask", dstate)
    ptr['wm_file'] = lo._get_aunique(dlayo, "wm_mask", dstate)
    # things that we want the names to create them later 
    ptr['dir_mask'] = lo._get_apth(dlayo, "sess_mask", dstate)
    ptr['mask_file'] = lo._get_apthglb(dlayo, "sess_mask", dstate)    
    #- Get runs' filenames and sort them 
    #-------------------------------------
    runs = [] 
    vals = range(1,nb_runs+1)
    for idx_run in vals:
        ds = dstate.copy()
        ds.update({dkeys["runs"]:idx_run})
        runs.append(sorted(lo._get_alistof(dlayo, "smoothed", ds)))
    ptr['runs'] = runs
    # Atlternative:
    #  runs = [ sorted(lo._get_alistof(dlayo, "smoothed", 
    #  lo.merge_two_dicts(dstate, {dkeys["runs"]:idx}))) for idx in range(1,nb_runs)]

    # compute_epi_mask(runs[0], opening=1, connected=True)
    #-----------------------------------------------------
    sess_mask = msk.compute_multi_epi_mask(runs, lower_cutoff=0.2, 
                    upper_cutoff=0.85, connected=True, opening=2, threshold=0.5)
    suf.rm_and_create(ptr['dir_mask'])
    if params['analysis']['apply_sess_mask']: sess_mask.to_filename(ptr['mask_file'])
    # TODO: check mask is reasonable - how ???

    # create the directory to write the extracted signals in
    save_is_true = params['analysis']['write_signals']
    if save_is_true: 
        suf.rm_and_create(ptr['signals_dir'])

    # - parameter file for condition names
    conditions = lo._get_aunique(dlayo, "conditions", dstate)
    with open(conditions) as fparam:
         sess_param = json.load(fparam)

    runs_info = {}
    for idx_run, run in enumerate(runs, 1): # /!\ starts at 1 /!\
        dstate[dkeys["runs"]] = idx_run
        ptr['file_names'] = run
        if verbose: print('\n' + '---'*9  + "\n" + "run{:02d}".format(idx_run))
        # TODO : fix this to have sess_param return motion['run1']='HIGH' etc
        ptr['motion'] = sess_param['motion'][idx_run-1] # sess_param['motion'] is 0 based

        runs_info["run{:02d}".format(idx_run)] = ptr #\
                    # do_one_run(ptr, dstate, dkeys, params, verbose=verbose)

    return runs_info

def do_one_run(run_curr, sess_curr, sub_curr, params, verbose=False):
    """
    """
    run_info = {}
    nvol = params['data']['nb_vol']
    dt = params['data']['TR']
    nb_run = params['data']['nb_run']

    file_names = run_curr['file_names']
    run_idx = run_curr['run_idx']
    run_idx0 = run_idx - 1
    assert run_idx0 >= 0
    assert run_idx0 < nb_run

    sess_idx = sess_curr['sess_idx']
    mvt_cond = run_curr['motion']
    dsig = sess_curr['dsig']
    mask = sess_curr['mask']

    low_freq = params['analysis']['filter']['low_freq']
    high_freq = params['analysis']['filter']['high_freq']
    
    # signal file names
    #-------------------
    _fn_sig = params['layout']['out']['signals']['signals+'] 
    # _fn_fsig = params['layout']['out']['signals']['f_signals+'] 
    fn_sig = osp.join(dsig, _fn_sig.format(sess_idx, run_idx) + mvt_cond)
    # fn_fsig = osp.join(dsig, _fn_fsig.format(run_idx)+mvt_cond)

    # extract signals and save them in preproc/roi_signals
    #-----------------------------------------------------
    min_vox_roi = params['analysis']['min_vox_in_roi']
    run_4d = concat_niimgs(file_names, ensure_ndim=4)
    signals, _issues, _info = ucr.extract_signals(sess_curr['droi'],
                                                  sess_curr['roi_prefix'], run_4d,
                                                  mask=mask, minvox=min_vox_roi)   
    # construct matrix of counfounds
    #-----------------------------------------------------
    arr_counf = []
    labs_counf = []
    #--- get WM 
    if params['analysis']['apply_wm']:
        wm_arr, wm_labs = ucr.extract_roi_run(
                            sess_curr['wm_dir'], sess_curr['wm_filename'], 
                            run_4d, check_lengh=nvol, verbose=verbose)
        labs_counf = labs_counf + wm_labs
        arr_counf.append(wm_arr)
        if verbose: print("applying wm \n")
    else: 
        wm_arr, wm_labs = None, None   
    #--- get CSF
    if params['analysis']['apply_csf']:
        csf_arr, csf_labs = ucr.extract_roi_run(
                            sess_curr['csf_dir'], sess_curr['csf_filename'], 
                            run_4d, check_lengh=nvol, verbose=verbose)
        labs_counf = labs_counf + csf_labs
        arr_counf.append(csf_arr)
        if verbose: print("applying csf \n")
    else: 
        csf_arr, csf_labs = None, None   
    #--- get GR 
    if params['analysis']['apply_global_reg']:
        gr_arr, gr_labs = ucr.extract_roi_run(
                            sess_curr['mask_dir'], sess_curr['mask_filename'], 
                            run_4d, check_lengh=nvol, verbose=verbose)
        labs_counf = labs_counf + gr_labs
        arr_counf.append(gr_arr)
        if verbose: print("applying GR \n")
    else: 
        gr_arr, gr_labs = None, None   
    #--- get MVT
    if params['analysis']['apply_mvt']:
        mvt_arr, mvt_labs = ucr.extract_mvt(sess_curr['mvtfile'], run_idx0, nvol, 
                                                                verbose=verbose)
        labs_counf = labs_counf + mvt_labs
        arr_counf.append(mvt_arr)
        if verbose: print("applying mvt \n")
    else: 
        mvt_arr, mvt_labs = None, None
    #--- get cosine functions;
    if params['analysis']['apply_filter']:
        bf_arr, bf_labs = ucr.extract_bf(low_freq, high_freq, nvol, dt, 
                                                                verbose=verbose)
        labs_counf = labs_counf + bf_labs
        arr_counf.append(bf_arr)
        if verbose: print("applying filter \n")
    else:
        bf_arr, bf_labs = None, None
    
    #--- put it together  
    # arr_counf = np.hstack((wm_arr, csf_arr, mvt_arr, bf_arr))
    # labs_counf = wm_labs + csf_labs + mvt_labs + bf_labs
    if arr_counf: 
        some_counfounds = True
        arr_counf = np.hstack(tuple(arr_counf))
    
    if verbose:
       #print("wm.shape {}, csf.shape {}, mvt.shape {}, bf.shape {}".format(
       #              wm_arr.shape, csf_arr.shape, mvt_arr.shape, bf_arr.shape))
       print(arr_counf.shape)
       print(labs_counf[:17])
    #run_info['shapes'] = (wm_arr.shape, csf_arr.shape, mvt_arr.shape, bf_arr.shape)
    #run_info['mean_csf'] = csf_arr.mean(axis=0)
    #run_info['mean_mvt'] = mvt_arr.mean(axis=0)

    # filter and save 
    #-----------------------------------------------------
    arr_sig, labels_sig = ucr._dict_signals_to_arr(signals) 
    if some_counfounds:
        arr_sig_f = ucr.R_proj(arr_counf, arr_sig)
    else:
        arr_sig_f = arr_sig

    run_info = dict(arr_sig=arr_sig, labels_sig=labels_sig, issues=_issues, 
                    info=_info, arr_sig_f=arr_sig_f, arr_counf=arr_counf, 
                    labs_counf=labs_counf)

    # save filtered signals 
    save_is_true = params['analysis']['write_signals']
    if save_is_true:
        np.savez(fn_sig, arr_sig=arr_sig, labels_sig=labels_sig, 
                         issues=_issues, info=_info, arr_sig_f=arr_sig_f,
                         arr_counf=arr_counf, labs_counf=labs_counf, params=params) 
        #np.savez(fn_fsig, arr_sig_f=arr_sig_f, arr_counf=arr_counf, labs_counf=labs_counf)

    return run_info

#------------------------------------------------------------------------------
# Running from command line
#------------------------------------------------------------------------------

if __name__ == "__main__":
    #
    parser = argparse.ArgumentParser()

    # Positional required arguments
    help_base_dir = "The directory containing sub01/sess??/... and json files: \n" + \
                    "analysis_parameters.json  data_parameters.json  directory_layout.json"
    parser.add_argument("base_directory", help=help_base_dir)

    # Optional arguments
    parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
    args = parser.parse_args()
    base_dir = args.base_directory
    verbose = args.verbose
    
    print("launching analyses on :", base_dir)
    assert osp.isdir(base_dir), '{} not a directory'.format(base_dir)

    info = process_all(base_dir, params=None, verbose=verbose)
   
    if verbose:
        print("\n------------------- Debug info ------------------ \n")
        print(info)
        print("\n------------------- Debug info ------------------ \n")

