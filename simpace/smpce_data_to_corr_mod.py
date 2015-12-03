import shutil
import os 
import os.path as osp
import glob as gb
import json
import numpy as np
#from six import string_types
#
import nilearn
from nilearn import masking as msk
from nilearn._utils import concat_niimgs
# from nilearn.image.image import _compute_mean
#
import utils._utils as ucr
import utils.setup_filenames as suf 
import argparse

import pdb
from six import string_types
import SimpleITK as sitk
import nibabel as nib

DIRLAYOUT = 'directory_layout.json'
DATAPARAM = 'data_parameters.json'
ANALPARAM = 'nuisance_parameters_6.json'

def get_params(dbase, verbose=False):
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
    ddata = params['data']

    # loop over subjects
    subj_idx = range(1,ddata['nb_sub']+1) # starts at 1, hence + 1
    subj_dirs = [osp.join(dbase, (dlayo['dir']['sub+']).format(idx)) for idx in subj_idx]
    # check all subj_dirs exists
    for sub_dir in subj_dirs:
        assert osp.isdir(sub_dir), "sub_dir"

    subjs_info = {}
    sub_curr = {}
    for sub_idx, sub_dir in enumerate(subj_dirs, 1): # start idx at 1
        sub_curr['sub_idx'] = sub_idx
        sub_curr['sub_dir'] = sub_dir
        sub_str = (dlayo['dir']['sub+']).format(sub_idx)
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
    nb_sess = params['data']['nb_sess']
    dlayo = params['layout']
    sess_idx = range(1, nb_sess+1)
    sess_dirs = [osp.join(sub_dir, (dlayo['dir']['sess+']).format(idx)) for idx in sess_idx]

    sesss_info = {} 
    sess_curr = {}
    for sess_idx, sess_dir in enumerate(sess_dirs, 1): # start idx at 1
        sess_curr['sess_idx'] = sess_idx
        sess_curr['sess_dir'] = sess_dir
        sess_str = (dlayo['dir']['sess+']).format(sess_idx)
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
    sub_idx = sub_curr['sub_idx']
    nb_runs = params['data']['nb_run'] 
    assert nb_runs == 4 # 4debug

    dlayo = params['layout']

    runs_dir = osp.join(sess_dir, dlayo['dir']['runs']) # should be preproc
    sess_curr['dir_runs'] = runs_dir

    if not 'data_prefix' in params['analysis']:
        dir_smooth_imgs = osp.join(runs_dir, dlayo['dir']['smooth'])
    else:
        dir_smooth_imgs = osp.join(runs_dir, params['analysis']['data_prefix']['type'])
    sess_curr['dir_smooth_imgs'] = dir_smooth_imgs #this is never used

    sess_curr['droi'] = osp.join(runs_dir, dlayo['atlas']['dir']) # 'registered_files'
    sess_curr['roi_prefix'] = dlayo['atlas']['prepat']            # 'rraal_*.nii'     
    sess_curr['dsig'] = osp.join(runs_dir, dlayo['out']['signals']['dir']) 

    save_is_true = params['analysis']['write_signals']
    if save_is_true: 
        # rm existing and recreate signal directory
        suf.rm_and_create(sess_curr['dsig'])

    sess_curr['dreal'] = osp.join(runs_dir, dlayo['dir']['realign'])

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
    pat_imgs_files = dlayo['pat']['sub+sess+run+']+"*.nii*"
                                # requires idx for sub, sess and run
    runs_pat = [pat_imgs_files.format(sub_idx, sess_idx, run_idx) \
                                        for run_idx in range(1, nb_runs+1)]
                                # /!\  start idx at 1 requires nb_runs+1  /!\
    runs = [gb.glob(osp.join(dir_smooth_imgs, pat)) for pat in runs_pat]
    # /!\ATTENTION:/!\ must sort the files with filename, should sort in time
    for run in runs: run.sort()
    sess_curr['runs'] = runs
    
    # compute session wide mask
    #-----------------------------------------------------
    # compute_epi_mask(runs[0], opening=1, connected=True)
    dir_mask = osp.join(runs_dir, dlayo['out']['sess_mask']['dir'])
    sess_curr['mask_dir'] = dir_mask
    sess_curr['mask_filename'] = dlayo['out']['sess_mask']['roi_mask']

    sess_mask = None
    # TODO : separate compute mask and apply
    if params['analysis']['apply_sess_mask']:
        mask_fullfile = osp.join(sess_curr['mask_dir'], sess_curr['mask_filename'])
        if not osp.isfile(mask_fullfile):
            sess_mask = msk.compute_multi_epi_mask(runs, lower_cutoff=0.2, 
                        upper_cutoff=0.85, connected=True, opening=2, threshold=0.5)
            suf.rm_and_create(dir_mask)
            sess_mask.to_filename(mask_fullfile)
        else:
            sess_mask = nib.load(mask_fullfile) #CORRECT? nib object?
    sess_curr['mask'] = sess_mask

    # TODO
    # check mask is reasonable - how ???

    # - mvt file
    # example : mvtfile = osp.join(dreal,'rp_asub01_sess01_run01-0006.txt')
    # /!\ will always be run01 for spm /!\
    mvtpat = dlayo['spm_mvt']['mvtrun1+'].format(sub_idx, sess_idx) 
    mvtfile = gb.glob(osp.join(sess_curr['dreal'], mvtpat))
    mvtfile = suf._check_glob_res(mvtfile, ensure=1, files_only=True)
    sess_curr['mvtfile'] = mvtfile

    # - parameter file for condition names
    param_pattern = (dlayo['pat']['sub+sess+']).format(sub_idx, sess_idx)
    paramfile = gb.glob(osp.join(sess_dir, param_pattern + "params*"))
    print sess_dir
    print param_pattern
    print osp.join(sess_dir, param_pattern + "params")
    print paramfile
    paramfile = suf._check_glob_res(paramfile, ensure=1, files_only=True)
    with open(paramfile) as fparam:
         sess_param = json.load(fparam)

    runs_info = {}
    run_curr = {}
    for idx_run, run in enumerate(runs, 1): # /!\ starts at 1 /!\
        run_curr['run_idx'] = idx_run
        run_curr['file_names'] = run
        if verbose: print('\n' + '---'*9  + "\n" + "run{:02d}".format(idx_run))
        # TODO : fix this to have sess_param return motion['run1']='HIGH' etc
        run_curr['motion'] = sess_param['motion'][idx_run-1] # sess_param['motion'] is 0 based

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

    file_names = run_curr['file_names']
    run_idx = run_curr['run_idx']
    run_idx0 = run_idx - 1
    assert run_idx0 >= 0
    assert run_idx0 < nb_run

    #sub_idx = sub_curr['sub_idx']
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

    run_4d = concat_niimgs(file_names, ensure_ndim=4)

    # construct matrix of counfounds
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


    # extract signals and save them in preproc/roi_signals
    #-----------------------------------------------------

    if params['analysis']['nuisance_level']['type'] == 'roi':
    
        min_vox_roi = params['analysis']['min_vox_in_roi']
        signals, _issues, _info = ucr.extract_signals(sess_curr['droi'], 
                                                      sess_curr['roi_prefix'], run_4d, 
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
        # save_is_true = params['analysis']['write_signals']
        # if save_is_true:
        #     np.savez(fn_sig, arr_sig=arr_sig, labels_sig=labels_sig, 
        #                      issues=_issues, info=_info, arr_sig_f=arr_sig_f,
        #                      arr_counf=arr_counf, labs_counf=labs_counf, params=params) 
            #np.savez(fn_fsig, arr_sig_f=arr_sig_f, arr_counf=arr_counf, labs_counf=labs_counf)

        return run_info


    elif params['analysis']['nuisance_level']['type'] == 'voxel':
        
        img_arr = np.asarray(run_4d.get_data())
        Nvols = img_arr.shape[3]
        Nvoxels = img_arr.size/Nvols
        img_arr_re = img_arr.reshape(Nvoxels,Nvols)
        img_arr_re = np.transpose(img_arr_re)

        if some_counfounds:
            arr_sig_f = ucr.R_proj(arr_counf, img_arr_re)
        else:
            arr_sig_f = img_arr_re
        
        arr_sig_f = np.transpose(arr_sig_f)

        arr_sig_f = arr_sig_f.reshape(img_arr.shape[0], img_arr.shape[1], img_arr.shape[2], Nvols)


        # if params['analysis']['write_signals']:

        save_prefix = params['analysis']['output_name']['type']

        save_dir = os.path.join( sess_curr['dir_runs'], save_prefix )
        
        if not os.path.isdir(save_dir):
            os.mkdir(save_dir)

        np.savetxt(os.path.join(save_dir, save_prefix + '_run_' + str(run_idx) + '.txt'), arr_counf)

        for ivol in range(0,Nvols):
            data = arr_sig_f[:,:,:,ivol]
            old_file = file_names[ivol]
            new_file = os.path.join(save_dir, save_prefix + '_' + os.path.basename(old_file) )
            
            old_obj = nib.load(old_file)
            affine = old_obj.get_affine()
            new_obj = nib.Nifti1Image(data, affine)
            nib.save(new_obj, new_file)

        affine = run_4d.get_affine()
        new_img = nib.Nifti1Image(arr_sig_f, affine)
        new_file = os.path.join(save_dir, save_prefix + '_run_' + str(run_idx) + '.nii.gz' )

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
    parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
    args = parser.parse_args()
    base_dir = args.base_directory
    verbose = args.verbose
    
    print("launching analyses on :", base_dir)
    print("parameters in ", ANALPARAM)
    assert osp.isdir(base_dir), '{} not a directory'.format(base_dir)

    info = process_all(base_dir, params=None, verbose=verbose)
   
    if verbose:
        print("\n------------------- Debug info ------------------ \n")
        print(info)
        print("\n------------------- Debug info ------------------ \n")

