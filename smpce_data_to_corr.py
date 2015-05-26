import shutil
import os 
import os.path as osp
import glob as gb
import json
import numpy as np
from six import string_types
#
from nilearn import masking as msk
from nilearn._utils import concat_niimgs
#
import utils._utils as ucr
import utils.setup_filenames as suf 

SKIPVOL = 4
RUNNBVOL = 200
VOLNB = RUNNBVOL - SKIPVOL
SUBNB = 1
SUBNBstr = '{:02d}'.format(SUBNB) 
TR = 2.0

DIRLAYOUT = 'directory_layout.json'
DATAPARAM = 'data_parameters.json'
ANALPARAM = 'analysis_parameters.json'

# from nilearn.image.image import _compute_mean


def file_or_dir(dlo, basedir, action):
    """
    dlo: directory layout and patterns, specific 
    """
    pass



def process_all(dbase):
    """
    parameters:
    -----------
    dbase:  string
            base directory containing subjects directory
            and jason files
    """

    # Read json files at the base directory
    fn_layo = osp.join(dbase, DIRLAYOUT)
    fn_data = osp.join(dbase, DATAPARAM)
    fn_anal = osp.join(dbase, ANALPARAM)

    with open(fn_layo) as flayo:
        dlayo = json.load(flayo)
    with open(fn_data) as fdata:
        ddata = json.load(fdata)
    with open(fn_anal) as fanal:
        danal = json.load(fanal)

    params = {'layout': dlayo, 'data_param': ddata, 'analy_param': danal}

    # loop over subjects
    subj_idx = range(1,ddata['nb_sub']+1) # starts at 1, hence + 1
    subj_dirs = [osp.join(dbase, (dlayo['dir']['sub+']).format(idx)) for idx in subj_idx]

    subjs_info = {}
    sub_curr = {}
    for sub_idx, sub_dir in enumerate(subj_dirs, 1): # start idx at 1
        sub_curr['sub_idx'] = sub_idx
        sub_curr['sub_dir'] = sub_dir
        subjs_info[sub_dir] =  do_one_subject(sub_curr, params) 

    return subjs_info

def do_one_subject(sub_curr, params):
    """
    Take the bdirectory at the moment subject
    and launch sessions processing

    parameters:
    -----------
    sub_dir: string
            subject base directory
    params: dict
            parameters for layout, data and analysis
            
    """
    sub_idx, sub_dir = sub_curr['sub_idx'], sub_curr['sub_dir']
    nb_sess = params['data_param']['nb_sess']
    dlayo = params['layout']
    sess_idx = range(1, nb_sess+1)
    sess_dirs = [osp.join(sub_dir, (dlayo['dir']['sess+']).format(idx)) for idx in sess_idx]

    sesss_info = {} 
    sess_curr = {}
    for sess_idx, sess_dir in enumerate(sess_dirs, 1): # start idx at 1
        sess_curr['sess_idx'] = sess_idx
        sess_curr['sess_dir'] = sess_dir
        sesss_info[sess_dir] = do_one_sess(sess_curr, sub_curr, params) 

    return sesss_info
    
def do_one_sess(sess_curr, sub_curr, params):
    """
    """

    sess_idx = sess_curr['sess_idx']
    sess_dir = sess_curr['sess_dir']
    sub_idx = sub_curr['sub_idx']
    nb_runs = params['data_param']['nb_run'] 
    assert nb_runs == 4 # 4debug

    dlayo = params['layout']

    runs_dir = osp.join(sess_dir, dlayo['dir']['runs']) # should be preproc
    sess_curr['dir_runs'] = runs_dir

    dir_smooth_imgs = osp.join(runs_dir, dlayo['dir']['smooth'])
    sess_curr['dir_smooth_imgs'] = dir_smooth_imgs

    droi = osp.join(runs_dir, dlayo['dir']['roi']) # 'registered_files'
    sess_curr['droi'] = droi
    sess_curr['roi_prefix'] = 'rraal_*.nii' # TODO: should be at project level
    
    dsig = osp.join(runs_dir, dlayo['dir']['signals']) #'extracted_signals'
    sess_curr['dsig'] = dsig
    # rm existing and recreate signal directory
    suf.rm_and_create(dsig)

    dreal = osp.join(runs_dir, dlayo['dir']['realign'])
    sess_curr['dreal'] = dreal

    sess_curr['csf_dir'] = osp.join(runs_dir, 'csf_mask')
    sess_curr['csf_filename'] =  'csf_map_final.nii.gz'

    #- Get runs' filenames
    #------------------------
    pat_imgs_files = dlayo['pat']['sub+sess+run+'] 
                                # will require idx for sub, sess and run
    runs_pat = [pat_imgs_files.format(sub_idx, sess_idx, run_idx) 
                                        for run_idx in range(1, nb_runs+1)]
                                # careful : start idx at 1 requires nb_runs+1
    runs = [gb.glob(osp.join(dir_smooth_imgs, pat)) for pat in runs_pat]
    for run in runs: run.sort()
    # ATTENTION: must sort the files - assume filenames will sort in time
    sess_curr['runs'] = runs
    
    # compute session wide mask
    #-----------------------------------------------------
    sess_mask = msk.compute_multi_epi_mask(runs, lower_cutoff=0.2, 
                    upper_cutoff=0.85, connected=True, opening=2, threshold=0.5)
    sess_curr['mask'] = sess_mask
    # TODO
    # check mask is reasonable - how ???
    # store sess mask 
    # compute_epi_mask(runs[0], opening=1, connected=True)

    # - mvt file
    # example : mvtfile = osp.join(dreal,'rp_asub01_sess01_run01-0006.txt')
    # will always be run01 for spm ?
    mvtpat = ('rp_asub{:02d}_sess{:02d}_run01' + '*-00*.txt').format(sub_idx, sess_idx) 
    mvtfile = gb.glob(osp.join(dreal, mvtpat))
    mvtfile = suf._check_glob_res(mvtfile, ensure=1, files_only=True)
    sess_curr['mvtfile'] = mvtfile

    # - parameter file for condition names
    param_pattern = (dlayo['pat']['sub+sess+']).format(sub_idx, sess_idx)
    paramfile = gb.glob(osp.join(sess_dir, param_pattern + "params"))
    paramfile = suf._check_glob_res(paramfile, ensure=1, files_only=True)
    with open(paramfile) as fparam:
         sess_param = json.load(fparam)

    runs_info = {}
    run_curr = {}
    for idx_run, run in enumerate(runs, 1): # start at 1;w
        run_curr['run_idx'] = idx_run
        run_curr['file_names'] = run
        run_curr['motion'] = sess_param['motion'][idx_run]

        runs_info["run_{:02d}".format(idx_run)] = \
                    do_one_run(run_curr, sess_curr, sub_curr, params)

    return runs_info

#def do_one_run(file_names, idx_run, mask, verbose=1):
def do_one_run(run_curr, sess_curr, sub_curr, params, verbose=1):
    """
    """
    run_info = {}
    nvol = params['data_param']['nb_vol']
    dt = params['data_param']['TR']

    file_names = run_curr['file_names']
    run_idx = run_curr['run_idx']
    #sub_idx = sub_curr['sub_idx']
    #sess_idx = sess_curr['sess_idx']
    mvt_cond = run_curr['motion']
    dsig = sess_curr['dsig']
    mask = sess_curr['mask']

    low_freq = params['analy_param']['low_freq']
    high_freq = params['analy_param']['high_freq']
    
    # file names
    #---------------
    # - signal files
    fn_sig = osp.join(dsig, 'signal_run{:02d}_'.format(run_idx) + mvt_cond)
    fn_fsig = osp.join(dsig, 'filtered_signal_run{:02d}_'.format(run_idx)+mvt_cond)

    # extract signals and save them in preproc/roi_signals
    #-----------------------------------------------------
    run_4d = concat_niimgs(file_names, ensure_ndim=4)
    signals, _issues, _info = ucr.extract_signals(sess_curr['droi'], 
                                                  sess_curr['roi_prefix'], run_4d, 
                                                  mask=mask, minvox=1)   
    arr_sig, labels_sig = ucr._dict_signals_to_arr(signals)
    np.savez(fn_sig, arr_sig, labels_sig) 

    # construct matrix of counfounds
    #-----------------------------------------------------
    #--- get CSF
    csf_arr, csf_labs = ucr.extract_roi_run(
                            sess_curr['csf_dir'], sess_curr['csf_filename'], 
                            run_4d, check_lengh=nvol, verbose=verbose)
    #--- get MVT
    mvt_arr, mvt_labs = ucr.extract_mvt(sess_curr['mvtfile'], run_idx, nvol, 
                                                                verbose=verbose)
    #--- get cosine functions;
    bf_arr, bf_labs = ucr.extract_bf(low_freq, high_freq, nvol, dt, 
                                                                verbose=verbose)
    #--- put it together  
    arr_counf = np.hstack((csf_arr, mvt_arr, bf_arr))
    labs_counf = csf_labs + mvt_labs + bf_labs
    if verbose:
       print("csf.shape {}, mvt.shape {}, bf.shape {}".format(
                     csf_arr.shape, mvt_arr.shape, bf_arr.shape))
    run_info['shapes'] = (csf_arr.shape, mvt_arr.shape, bf_arr.shape)

    # filter and compute correlation
    #-----------------------------------------------------
    arr_sig_f = ucr.R_proj(arr_counf, arr_sig)

    # save filtered signals 
    np.savez(fn_fsig, arr_sig_f, labels_sig, arr_counf, labs_counf)

    return run_info


