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

# from nilearn.image.image import _compute_mean


def file_or_dir(dlo, basedir, action):
    """
    dlo: directory layout and patterns, specific 
    """
    pass



def do_one_sess(DDIR, idx_sess):
    """
    take a session directory and extract the roi signals
    """
    runs_pat = ['*sess{:02d}_run{:02d}-0*'.format(idx_sess, idx_run) 
                                                for idx_run in [1,2,3,4]]
    runs = [gb.glob(osp.join(DSRA, pat)) for pat in runs_pat]
    for run in runs: run.sort()

    # compute session wide mask
    #-----------------------------------------------------

    sess_mask = msk.compute_multi_epi_mask(runs, lower_cutoff=0.2, 
                    upper_cutoff=0.85, connected=True, opening=2, threshold=0.5)

    #compute_epi_mask(runs[0], opening=1, connected=True)
    
    # check mask is reasonable - how ???

    for idx_run, run in enumerate(runs):
        info = do_one_run(run, idx_run, sess_mask)


def do_one_run(file_names, idx_run, mask, verbose=1):
    """
    """
    #substr = SUBNBstr
    #sesstr = 'sess{:02d}'.format(idx_sess)
    #runstr = 'run{:02d}'.format(idx_run)
    # extract base directory from first file
    ddir = osp.dirname(osp.dirname(file_names[0])) # preproc dir
    dsess = osp.dirname(ddir) # session directory is one up
    droi = osp.join(ddir, 'registered_files')
    dsig = osp.join(ddir, 'extracted_signals')
    dreal = osp.join(ddir, 'realign')
    # - parameter (condition name)
    paramfile = gb.glob(osp.join(dsess,"sub??_sess??_params"))
    paramfile = suf._check_glob_res(paramfile, ensure=1, files_only=True)
    with open(paramfile) as fparam:
         param = json.load(fparam)
         mvt_cond = str(param['motion'][idx_run])
    # file names
    #---------------
    # - signal files
    fn_sig = osp.join(dsig, 'signal_run{:02d}_'.format(idx_run) + mvt_cond)
    fn_fsig = osp.join(dsig, 'filtered_signal_run{:02d}_'.format(idx_run)+mvt_cond)
    # - csf file
    dcsf = osp.join(ddir, 'csf_mask')
    csf_filename =  'csf_map_final.nii.gz'
    # - mvt file
    # example : mvtfile = osp.join(dreal,'rp_asub01_sess01_run01-0006.txt')
    mvtpat = ('rp_asub??_sess??_run{:02d}' + '*-00*.txt').format(idx_run+1)
    mvtfile = gb.glob(osp.join(dreal, mvtpat))
    mvtfile = suf._check_glob_res(mvtfile, ensure=1, files_only=True)
    # - atlas roi files
    roi_prefix = 'rraal_*.nii'
    #- other parameters:
    #--------------------
    low_freq = 0.01
    high_freq = 0.1
    dt = TR

    # extract signals and save them in preproc/roi_signals
    #-----------------------------------------------------
    run_4d = concat_niimgs(file_names, ensure_ndim=4)
    signals, _issues, _info = ucr.extract_signals(droi, roi_prefix, run_4d, 
                                                        mask=mask, minvox=1)   
    arr_sig, labels_sig = ucr._dict_signals_to_arr(signals)

    # rm previous directory
    suf.rm_and_create(dsig)
    np.savez(fn_sig, arr_sig, labels_sig) 

    # construct matrix of counfounds
    #-----------------------------------------------------
    #--- get CSF
    csf_arr, csf_labs = ucr.extract_roi_run(dcsf, csf_filename, run_4d, 
                                             check_lengh=VOLNB, verbose=verbose)
    #--- get MVT
    mvt_arr, mvt_labs = ucr.extract_mvt(mvtfile, idx_run, VOLNB, verbose=verbose)
    #--- get cosine functions;
    bf_arr, bf_labs = ucr.extract_bf(low_freq, high_freq, VOLNB, dt, 
                                                            verbose=verbose)
    #--- put it together  
    arr_counf = np.hstack((csf_arr, mvt_arr, bf_arr))
    labs_counf = csf_labs + mvt_labs + bf_labs
    if verbose:
       print("csf.shape {}, mvt.shape {}, bf.shape {}".format(
                     csf_arr.shape, mvt_arr.shape, bf_arr.shape))

    # filter and compute correlation
    #-----------------------------------------------------
    arr_sig_f = ucr.R_proj(arr_counf, arr_sig)

    # save filtered signals 
    np.savez(fn_fsig, arr_sig_f, labels_sig, arr_counf, labs_counf) 


