import shutil
import os 
import os.path as osp
import glob as gb
import json
import numpy as np
import math
#from six import string_types
#
from nilearn import masking as msk
from nilearn._utils import concat_niimgs
# from nilearn.image.image import _compute_mean

import utils._utils as ucr
import utils.setup_filenames as suf 
import argparse

import nibabel as nib

DIRLAYOUT = 'directory_layout_bids.json'
DATAPARAM = 'data_parameters.json'
# ANALPARAM = 'nuisance_parameters_0.json'
# ANALPARAM = 'analysis_parameters_extract.json'


def get_params(dbase, trajectory, sessions, ANALPARAM, verbose=False):
    """
    parameters:
    -----------
    dbase:  string
            base directory containing subjects directory
            and jason files
    """
    json_dir = '/home/despo/arielle/simpace/simpace'

    # Read json files at the base directory
    fn_layo = osp.join(json_dir, DIRLAYOUT)
    fn_data = osp.join(json_dir, DATAPARAM)
    fn_anal = osp.join(json_dir, ANALPARAM)

    with open(fn_layo) as flayo:
        dlayo = json.load(flayo)
    with open(fn_data) as fdata:
        ddata = json.load(fdata)
    with open(fn_anal) as fanal:
        danal = json.load(fanal)

    if not sessions:
        ddata['sessions'] = '*'
    else:
        ddata['sessions'] = sessions
    ddata['trajectory'] = trajectory

    params = {'layout': dlayo, 'data': ddata, 'analysis': danal}
    return params


def process_all(dbase, subject, trajectory, sessions, ANALPARAM, params=None, verbose=False):
    """
    parameters:
    -----------
    dbase:  string
            base directory containing subjects directory
            and jason files
    """
    if not params:
        params = get_params(dbase, trajectory, sessions, ANALPARAM, verbose=verbose)
    
    dlayo = params['layout']
    ddata = params['data']

    # set up subject (not loop over subjects)
    # subj_idx = subject  # input to script
    subj_dir = osp.join(dbase, (dlayo['dir']['sub+']).format(subject))
    # check subj_dir exists
    assert osp.isdir(subj_dir), "sub_dir"

    sub_curr = {}
    subjs_info = {}
    sub_curr['sub_idx'] = subject  # input to script
    sub_curr['sub_dir'] = subj_dir
    sub_str = (dlayo['dir']['sub+']).format(subject)

    subjs_info[sub_str] =  do_one_subject(sub_curr, params, verbose=verbose)

    return subjs_info


def do_one_subject(sub_curr, params, verbose=False):
    """
    launch sessions processing for sub_curr 

    parameters:
    -----------
    sub_curr: dict 
            contains subject base directory
            contains subject index 
    params: dict
            parameters for layout, data and analysis
            
    """
    sub_idx, sub_dir = sub_curr['sub_idx'], sub_curr['sub_dir']
    # nb_sess = params['data']['nb_sess']
    dlayo = params['layout']
    traj = params['data']['trajectory']
    sess_list = params['data']['sessions']
    if sess_list is '*':
        sess_dirs = gb.glob(osp.join(sub_dir, '*' + traj + '*'))
    else:
        sess_dirs = [osp.join(sub_dir, (dlayo['dir']['sess+']).format(traj, idx)) for idx in sess_list]

    sesss_info = {} 
    sess_curr = {}
    for sess_idx, sess_dir in enumerate(sess_dirs, 1): # start idx at 1
        curr_sess = sess_list[sess_idx-1]
        sess_curr['sess_idx'] = curr_sess
        sess_curr['sess_dir'] = sess_dir
        sess_curr['traj'] = traj
        sess_str = (dlayo['dir']['sess+']).format(traj, curr_sess)
        if verbose: print('\n' + '---'*11  + "\n" + sess_str)
        sesss_info[sess_str] = do_one_sess(sess_curr, sub_curr, params, verbose=verbose) 

    return sesss_info


def do_one_sess(sess_curr, sub_curr, params, verbose=False):
    """
    launch runs processing for sess_curr 

    parameters:
    -----------
    sess_curr: dict 
            contains sess base directory
            contains sess index 
    params: dict
            parameters for layout, data and analysis
    """            

    sess_idx = sess_curr['sess_idx']
    sess_dir = sess_curr['sess_dir']
    traj = sess_curr['traj']
    sub_idx = sub_curr['sub_idx']
    nb_runs = params['data']['nb_run'] 
    assert nb_runs == 4 # 4debug

    dlayo = params['layout']

    runs_dir = osp.join(sess_dir, dlayo['dir']['runs']) # preproc
    sess_curr['dir_runs'] = runs_dir

    # specify directory of data to be extracted from - smooth is default as specified in directory_layout.json
    if 'data_prefix' in params['analysis']:
        params['layout']['dir']['smooth'] = params['analysis']['data_prefix']['dir']

    dir_smooth_imgs = osp.join(runs_dir, params['layout']['dir']['smooth'])

    sess_curr['droi'] = osp.join(runs_dir, dlayo['atlas']['dir'])  # e.g. preproc/coreg
    sess_curr['roi_name'] = dlayo['atlas']['name']  # names of atlas(e)s
    sess_curr['roi_prefix'] = dlayo['atlas']['prepat']  # '*.nii'
    sess_curr['dsig'] = osp.join(runs_dir, dlayo['out']['signals']['dir'])  # extracted_signals dir

    save_is_true = params['analysis']['write_signals']

    sess_curr['dreal'] = osp.join(runs_dir, dlayo['dir']['realign'])  # used?

    #- csf dir and file
    sess_curr['csf_dir'] = osp.join(runs_dir, dlayo['csf']['dir'])
    csf_file = gb.glob(osp.join(sess_curr['csf_dir'], dlayo['csf']['roi_mask']))
    if not csf_file: print("glob empty: {} {}".format(
                                sess_curr['csf_dir'], dlayo['csf']['roi_mask']))
    csf_file = suf._check_glob_res(csf_file, ensure=1, files_only=True)
    sess_curr['csf_filename'] =  dlayo['csf']['roi_mask']

    #- wm dir and file
    sess_curr['wm_dir'] = osp.join(runs_dir, dlayo['wm']['dir'])
    wm_file = gb.glob(osp.join(sess_curr['wm_dir'], dlayo['wm']['roi_mask']))
    if not wm_file: print("glob empty: {} {}".format(
                                sess_curr['wm_dir'], dlayo['wm']['roi_mask']))
    wm_file = suf._check_glob_res(wm_file, ensure=1, files_only=True)
    sess_curr['wm_filename'] =  dlayo['wm']['roi_mask']

    #- Get runs' filenames
    #------------------------
    pat_files = (dlayo['pat']['sub+sess+run+']).format(sub_idx, traj, sess_idx) + "*.nii*"
    runs = gb.glob(osp.join(dir_smooth_imgs, pat_files))
    print runs
    print dir_smooth_imgs
    print pat_files

    # If specified pattern doesn't exist, used base prefix specified in analysis parameters json file
    # *Need to update w/ BIDS formatting to use this functionality
    if not runs:
        pat_files = '*' + params['analysis']['data_prefix']['base'] + '*.nii*'
        runs = gb.glob(osp.join(dir_smooth_imgs, pat_files))

        # for irun in range(1, nb_runs+1):
        #     runs_pat[irun-1] = '*' + params['analysis']['data_prefix']['base'] + '*_run0' + str(irun) + '*.nii*'

    # /!\ATTENTION:/!\ sorting not done w/ BIDS formatting, glob default sorts in time
    sess_curr['runs'] = runs
    print runs
    
    # compute session wide mask
    #-----------------------------------------------------
    # compute_epi_mask(runs[0], opening=1, connected=True)
    dir_mask = osp.join(runs_dir, dlayo['out']['sess_mask']['dir'])
    sess_curr['mask_dir'] = dir_mask
    sess_curr['mask_filename'] = dlayo['out']['sess_mask']['roi_mask']

    sess_mask = None  # don't think that this does anything
    # TODO : separate compute mask and apply
    if params['analysis']['apply_sess_mask']:
        mask_fullfile = osp.join(sess_curr['mask_dir'], sess_curr['mask_filename'])
        if not osp.isfile(mask_fullfile):
            sess_mask = msk.compute_multi_epi_mask(runs, lower_cutoff=0.2, 
                        upper_cutoff=0.85, connected=True, opening=2, threshold=0.5)
            suf.rm_and_create(dir_mask)
            sess_mask.to_filename(mask_fullfile)
        else:
            sess_mask = nib.load(mask_fullfile)

    sess_curr['mask'] = sess_mask

    # TODO
    # check mask is reasonable - how ???

    runs_info = {}
    run_curr = {}
    for idx_run, run in enumerate(runs, 1): # /!\ starts at 1 /!\
        run_curr['run_idx'] = idx_run
        run_curr['file_names'] = run
        run_base = os.path.basename(run)

        if 'data_prefix' in params['analysis']:
            idx_st = run_base.index(params['analysis']['data_prefix']['dir']) + \
                len(params['analysis']['data_prefix']['dir'])+1
            idx_end = run_base.index('.nii')
        else:
            idx_st = run_base.index(dlayo['pat']['run+']) + len(dlayo['pat']['run+'])
            idx_end = run_base.index(dlayo['pat']['postfix'])

        if verbose: print('\n' + '---'*9  + "\n" + "run{:02d}".format(idx_run))
        # TODO : fix this to have sess_param return motion['run1']='HIGH' etc
        run_curr['motion'] = run_base[idx_st:idx_end]

        mvtpat = dlayo['spm_mvt']['mvtrun'].format(sub_idx, traj, sess_idx, run_curr['motion'])
        mvtfile = gb.glob(osp.join(sess_curr['dreal'], mvtpat))

        if not mvtfile:
            mvtfile = ''
        else:
            mvtfile = suf._check_glob_res(mvtfile, ensure=1, files_only=True)
        run_curr['mvtfile'] = mvtfile

        runs_info["run{:02d}".format(idx_run)] = \
                    do_one_run(run_curr, sess_curr, sub_curr, params, verbose=verbose)

    return runs_info


def do_one_run(run_curr, sess_curr, sub_curr, params, verbose=False):
    """
    """
    run_info = {}
    nvol = params['data']['nb_vol']
    dt = params['data']['TR']
    nb_run = params['data']['nb_run']

    file_name = run_curr['file_names']
    run_idx = run_curr['run_idx']
    run_idx0 = run_idx - 1
    assert run_idx0 >= 0
    assert run_idx0 < nb_run

    sub_idx = sub_curr['sub_idx']
    sess_idx = sess_curr['sess_idx']
    traj = sess_curr['traj']
    mvt_cond = run_curr['motion']
    mvt_file = run_curr['mvtfile']
    dsig = sess_curr['dsig']
    mask = sess_curr['mask']

    save_is_true = params['analysis']['write_signals']

    # Load in data - assumes 4D format
    run_4d = nib.load(file_name)

    # construct matrix of confounds
    #-----------------------------------------------------
    arr_counf = []
    labs_counf = []

    if params['analysis']['apply_wm'] or params['analysis']['apply_csf']:
        if params['analysis']['wm']['type'] == 'PC' or params['analysis']['csf']['type'] == 'PC':
            n_PC = int(params['analysis']['PC_number']['type'])
        else:
            n_PC = 5

    #--- get WM 
    if params['analysis']['apply_wm']:
        wm_arr, wm_labs = ucr.extract_roi_run(
                            sess_curr['wm_dir'], sess_curr['wm_filename'], 
                            run_4d, check_lengh=nvol, verbose=verbose,
                            signal_type=params['analysis']['wm']['type'], 
                            n_PC=n_PC)
        labs_counf = labs_counf + wm_labs
        arr_counf.append(wm_arr)
        if verbose: print("applying wm \n")
    else: 
        wm_arr, wm_labs = None, None

    #--- get CSF
    if params['analysis']['apply_csf']:
        csf_arr, csf_labs = ucr.extract_roi_run(
                            sess_curr['csf_dir'], sess_curr['csf_filename'], 
                            run_4d, check_lengh=nvol, verbose=verbose,
                            signal_type=params['analysis']['csf']['type'], 
                            n_PC=n_PC)
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
        mvt_arr, mvt_labs = ucr.extract_mvt_perrun(mvt_file, nvol,
                                                                verbose=verbose)
        labs_counf = labs_counf + mvt_labs
        arr_counf.append(mvt_arr)
        if verbose: print("applying mvt \n")
    else: 
        mvt_arr, mvt_labs = None, None

    #--- get cosine functions;
    if params['analysis']['apply_filter']:
        low_freq = float(params['analysis']['filter']['low_freq'])
        high_freq = float(params['analysis']['filter']['high_freq'])

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
    else:
        some_counfounds = False
        print 'Not regressing anything out of data!'
    
    if verbose:
       #print("wm.shape {}, csf.shape {}, mvt.shape {}, bf.shape {}".format(
       #              wm_arr.shape, csf_arr.shape, mvt_arr.shape, bf_arr.shape))
       print(arr_counf.shape)
       print(labs_counf[:17])
    #run_info['shapes'] = (wm_arr.shape, csf_arr.shape, mvt_arr.shape, bf_arr.shape)
    #run_info['mean_csf'] = csf_arr.mean(axis=0)
    #run_info['mean_mvt'] = mvt_arr.mean(axis=0)


    # extract signals and save them in preproc/roi_signals
    #-----------------------------------------------------

    if params['analysis']['nuisance_level']['type'] == 'roi':

        atlas_names = (params['layout']['atlas']['name'])

        _fn_sig = params['layout']['out']['signals']['signals+']
        min_vox_roi = params['analysis']['min_vox_in_roi']

        if not osp.isdir(dsig):
            os.mkdir(dsig)

        for iatlas in atlas_names:

            # signal file names
            # _fn_fsig = params['layout']['out']['signals']['f_signals+']
            fn_sig = osp.join(dsig, params['layout']['dir']['smooth'] + "_" \
                + _fn_sig.format(sub_idx, traj, sess_idx, mvt_cond) + "_" + iatlas)
            # fn_fsig = osp.join(dsig, _fn_fsig.format(run_idx)+mvt_cond)

            print fn_sig

            roi_files = '*' + iatlas + sess_curr['roi_prefix']
            print roi_files
            signals, _issues, _info = ucr.extract_signals(sess_curr['droi'],
                                                          roi_files, run_4d,
                                                          mask=mask, minvox=min_vox_roi)
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
            if save_is_true:

                np.savez(fn_sig, arr_sig=arr_sig, labels_sig=labels_sig,
                                 issues=_issues, info=_info, arr_sig_f=arr_sig_f,
                                 arr_counf=arr_counf, labs_counf=labs_counf, params=params)
                # np.savez(fn_fsig, arr_sig_f=arr_sig_f, arr_counf=arr_counf, labs_counf=labs_counf)

        return run_info

    elif params['analysis']['nuisance_level']['type'] == 'voxel':
        
        img_arr = np.asarray(run_4d.get_data())
        Nvols = img_arr.shape[3]
        Nvoxels = img_arr.size/Nvols
        img_arr_re = img_arr.reshape(Nvoxels, Nvols)
        img_arr_re = np.transpose(img_arr_re)

        if some_counfounds:
            arr_sig_f = ucr.R_proj(arr_counf, img_arr_re)
        else:
            arr_sig_f = img_arr_re
            print 'This pipeline does not make sense!'
        
        arr_sig_f = np.transpose(arr_sig_f)

        arr_sig_f = arr_sig_f.reshape(img_arr.shape[0], img_arr.shape[1], img_arr.shape[2], Nvols)

        save_prefix = params['analysis']['output_name']['type']

        save_dir = os.path.join(sess_curr['dir_runs'], save_prefix)
        
        if not os.path.isdir(save_dir):
            os.mkdir(save_dir)

        save_name = os.path.join(save_dir, save_prefix + '_' + mvt_cond)
        np.savetxt(save_name + '.txt', arr_counf)

        #  write out a 4-d file
        affine = run_4d.get_affine()
        new_img = nib.Nifti1Image(arr_sig_f, affine)
        new_file = save_name + '.nii.gz'

        nib.save(new_img, new_file)


    

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
    parser.add_argument("subject_number")
    parser.add_argument("trajectory")
    parser.add_argument("session_number")
    parser.add_argument("analysis_params")
    parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
    args = parser.parse_args()
    base_dir = args.base_directory
    subject = int(args.subject_number)
    session = args.session_number
    if session is not '*':
        session = eval(session)
    trajectory = args.trajectory
    ANALPARAM = args.analysis_params
    verbose = args.verbose
    
    print("launching analyses on :", base_dir)
    print("parameters in ", ANALPARAM)
    assert osp.isdir(base_dir), '{} not a directory'.format(base_dir)

    info = process_all(base_dir, subject, trajectory, session, ANALPARAM, params=None, verbose=verbose)
   
    if verbose:
        print("\n------------------- Debug info ------------------ \n")
        print(info)
        print("\n------------------- Debug info ------------------ \n")

