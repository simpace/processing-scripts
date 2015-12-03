from __future__ import print_function
import os
import os.path as osp
import json
import numpy as np
import argparse
#from six import string_types
from nilearn import masking as msk
from nilearn._utils import concat_niimgs
import nibabel as nib
# from nilearn.image.image import _compute_mean
import utils._utils as ucr
import utils.setup_filenames as suf 
import layout as lo

DIRLAYOUT = 'directory_layout.json'
DATAPARAM = 'data_parameters.json'
ANALPARAM = 'analysis_parameters.json'


def get_params(dbase, addjson=None, verbose=False):
    """
    Read json files at the base directory dbase
    
    parameters:
    -----------
    dbase:  string
            base directory containing subjects directory
            *and* jason files
    addjson: string
            filename of an addtional json file
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

    if addjson:
        # make sure addjson is a json file 
        if not osp.isfile(addjson):
            addjson = osp.join(dbase, osp.basename(addjson)+'.json')
        if not osp.isfile(addjson):
            raise ValueError("json file does not exist : {} ".format(addjson))
        # read and populate params - assume 2 levels
        with open(addjson) as fjson:
            djson = json.load(fjson)
        for top_key in djson:
            for second_level_key in djson[top_key]:
                params[top_key][second_level_key] = djson[top_key][second_level_key]
                if verbose: print(top_key, second_level_key, djson[top_key][second_level_key])

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
    # do only those specified in analysis
    if params['analysis']['subjects']:
        _0idx = params['analysis']['subjects']
        idxs = [idxs[idx] for idx in _0idx]
        sub_dirs = [sub_dirs[idx] for idx in _0idx]

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
    if params['analysis']['sessions']:
        _0idx = params['analysis']['sessions']
        idxs = [idxs[idx] for idx in _0idx]
        sess_dirs = [sess_dirs[idx] for idx in _0idx]

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
    ptr['signals_dir'] = lo._get_apth(dlayo, "signals", dstate)
    #---------
    ptr['aal_dir'] = lo._get_apth(dlayo, "aal_roi", dstate) #  preproc/..
    ptr['aal_files'] = lo._get_alistof(dlayo, "aal_roi", dstate) #  
    ptr['aal_glb'] = dlayo['aal_roi']['glb']
    #---------
    ptr['gm_dir'] = lo._get_apth(dlayo, "gm_mask", dstate)
    ptr['gm_file'] = lo._get_aunique(dlayo, "gm_mask", dstate)
    ptr['gm_glb'] = dlayo['gm_mask']['glb']
    #--------- things that we want the names to create them later 
    ptr['mask_dir'] = lo._get_apth(dlayo, "sess_mask", dstate)
    # ptr['mask_file'] = lo._get_aunique(dlayo, "sess_mask", dstate)    
    ptr['mask_glb'] = dlayo['sess_mask']['glb']

    #- Get runs' filenames and sort them 
    #-------------------------------------
    runs = [] 
    run_indices = range(1,nb_runs+1)
    for idx_run in run_indices:
        ds = dstate.copy()
        ds.update({dkeys["runs"]:idx_run})
        runs.append(sorted(lo._get_alistof(dlayo, "smoothed", ds)))

    for a_run in runs:
        if not a_run:
            helpValueError = "\nPlease check what is in {}".format(
                                        lo._get_apth(dlayo, "smoothed", ds))
            helpValueError += print(lo._get_alistof(dlayo, "smoothed", ds))
            helpValueError += print(lo._get_apthglb(dlayo, "smoothed", ds, glbo=False))
            raise ValueError("runs are empty: {}".format(runs) + helpValueError)

    for a_run in runs:
        ucr.check_affines(a_run, verbose=verbose)

    ptr['runs'] = runs

    # Atlternative:
    #  runs = [ sorted(lo._get_alistof(dlayo, "smoothed", 
    #  lo.merge_two_dicts(dstate, {dkeys["runs"]:idx}))) for idx in range(1,nb_runs)]

    # compute_epi_mask(runs[0], opening=1, connected=True)
    #-----------------------------------------------------
    
    # the name of the file is always defined ... but not necessarily computed or written
    ptr['mask_file'] = osp.join(ptr['mask_dir'], ptr['mask_glb'])
    if params['analysis']['compute_sess_mask']:
        sess_mask = msk.compute_multi_epi_mask(runs, lower_cutoff=0.2, 
                    upper_cutoff=0.85, connected=True, opening=3, threshold=0.5)
        if params['analysis']['write_sess_mask']: 
            suf.rm_and_create(ptr['mask_dir'])
            sess_mask.to_filename(ptr['mask_file'])
            if verbose: print("checking mask affine: \n")
            ucr.check_affines([ptr["mask_file"]] + [r[0] for r in runs], verbose=verbose)

    # TODO: check mask is reasonable - how ???

    # create the directory to write the extracted signals in
    save_is_true = params['analysis']['write_signals']
    
    if save_is_true: 
        if not osp.exists(ptr['signals_dir']):
            os.makedirs(ptr['signals_dir'])
            if verbose: print('\n making directory {}'.format(ptr['signals_dir']))

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

        runs_info["run{:02d}".format(idx_run)] = ptr # save parameters 
        do_one_run(ptr, dstate, dkeys, params, verbose=verbose)

    return runs_info

def do_one_run(ptr,  dstate, dkeys, params, verbose=False):  # run_curr, sess_curr, sub_curr,
    """
    """
    run_info = {}
    dlayo = params['layout']
    nvol = params['data']['nb_vol']
    dt = params['data']['TR']
    nb_run = params['data']['nb_run']

    file_names = ptr['file_names']
    run_idx = dstate[dkeys["runs"]] #run_curr['run_idx']
    run_idx0 = run_idx - 1
    assert run_idx0 >= 0
    assert run_idx0 < nb_run

    # sess_idx = dstate[dkeys["sessions"]] # sess_curr['sess_idx']
    mvt_cond = ptr['motion'] # run_curr['motion']

    if params['analysis']['apply_sess_mask']:
        mask = ptr['mask_file'] # sess_curr['mask']
    else:
        mask = None

    low_freq = params['analysis']['filter']['low_freq']
    high_freq = params['analysis']['filter']['high_freq']
    
    # signal file names
    #-------------------
    fn_sig = lo._get_apthglb(dlayo, "signals", dstate, glob=False) + mvt_cond

    # extract signals and save them in preproc/roi_signals
    #-----------------------------------------------------
    min_vox_roi = params['analysis']['min_vox_in_roi']

    # if there is only one file per run, check it is a 4d file
    should_be_4d = False
    if len(file_names) == 1:
            img4d = nib.load(file_names[0])
            assert img4d.shape[-1] > 1, "img4d.shape = {}".format(img4d.shape)
            should_be_4d = True
    print("should_be_4d",should_be_4d, len(file_names))
    # concat if not already 4d
    if not should_be_4d:
        run_4d = concat_niimgs(file_names, ensure_ndim=4)
    else:
        run_4d = file_names[0]

    # Extract signals from ROI
    signals, _issues, _info = ucr.extract_signals(ptr['aal_dir'], ptr['aal_glb'],  
                                        run_4d, mask=mask, minvox=min_vox_roi, verbose=verbose)   
    # construct matrix of counfounds
    #-----------------------------------------------------
    arr_counf = []
    labs_counf = []
    #--- get WM 
    if params['analysis']['apply_wm']:
        wm_dir = lo._get_apth(dlayo, "wm_mask", dstate)
        #ptr['wm_file'] = lo._get_aunique(dlayo, "wm_mask", dstate)
        wm_glb = dlayo['wm_mask']['glb']
        #---------
        wm_arr, wm_labs = ucr.extract_roi_run( wm_dir, wm_glb, 
                            run_4d, check_lengh=nvol, verbose=verbose)
        labs_counf = labs_counf + wm_labs
        arr_counf.append(wm_arr)
        if verbose: print("applying wm \n")
    else: 
        wm_arr, wm_labs = None, None   
    #--- get CSF
    if params['analysis']['apply_csf']:
        csf_dir = lo._get_apth(dlayo, "csf_mask", dstate)
        #ptr['csf_file'] = lo._get_aunique(dlayo, "csf_mask", dstate)
        csf_glb = dlayo['csf_mask']['glb']
        #---------
        csf_arr, csf_labs = ucr.extract_roi_run(csf_dir, csf_glb, 
                                run_4d, check_lengh=nvol, verbose=verbose)
        labs_counf = labs_counf + csf_labs
        arr_counf.append(csf_arr)
        if verbose: print("applying csf \n")
    else: 
        csf_arr, csf_labs = None, None   
    #--- get GR 
    if params['analysis']['apply_global_reg']:
        gr_arr, gr_labs = ucr.extract_roi_run(
                            ptr['mask_dir'], ptr['mask_glb'], 
                            run_4d, check_lengh=nvol, verbose=verbose)
        labs_counf = labs_counf + gr_labs
        arr_counf.append(gr_arr)
        if verbose: print("applying GR \n")
    else: 
        gr_arr, gr_labs = None, None   
    #--- get MVT
    if params['analysis']['apply_mvt']:
        mvt_file = lo._get_aunique(dlayo, "mvt6params", dstate)
        mvt_arr, mvt_labs = ucr.extract_mvt(mvt_file, run_idx0, nvol, 
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
    some_counfounds = False
    if arr_counf: 
        some_counfounds = True
        arr_counf = np.hstack(tuple(arr_counf))
    
    if verbose:
       #print("wm.shape {}, csf.shape {}, mvt.shape {}, bf.shape {}".format(
       #              wm_arr.shape, csf_arr.shape, mvt_arr.shape, bf_arr.shape))
       if some_counfounds:
           print(arr_counf.shape)
           print(labs_counf[:np.min([17,arr_counf.shape[1]])])
       else: 
           print('no counfounds')

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
        print("\saving {}".format(fn_sig))
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

    tmpfile = osp.join(base_dir, 'tmp.log')
    with  open(tmpfile, 'w') as file_handle:
        print("\n------------------- Debug info ------------------ \n")
        print(info, file=file_handle)
        print("\n------------------- Debug info ------------------ \n")





