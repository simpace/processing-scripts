import shutil
import os 
import os.path as osp
import glob as gb
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

    sess_mask = msk.compute_multi_epi_mask(runs, lower_cutoff=0.2, upper_cutoff=0.85,
                           connected=True, opening=2, threshold=0.5)
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
    ddir = osp.dirname(file_names[0]) # preproc dir
    dsess = osp.dirname(ddir) # session directory is one up
    droi = osp.join(ddir, 'registered_files')
    dsig = osp.join(ddir, 'extracted_signals')
    fn_sig = osp.join(dsig, 'signal_run_{:02d}'.format(idx_run))
    fn_fsig = osp.join(dsig, 'flt_signal_run_{:02d}'.format(idx_run))
    # - csf file
    dcsf = osp.join(ddir, 'csf_mask')
    csf_filename =  'csf_map_final.nii.gz'
    # - mvt file
    dreal = osp.join(ddir, 'realign')
    mvtpat = 'rp_asub??_' + '*-0*.txt'
    #example : mvtfile = osp.join(dreal,'rp_asub01_sess01_run01-0006.txt')
    mvtfile = gb.glob(osp.join(dreal, mvtpat))
    mvtfile = suf._check_glob_res(mvtfile, ensure=1, files_only=True)
    # - atlas roi files
    roi_prefix = 'rraal_*.nii'

    #- other parameters:
    low_freq = 0.01
    high_freq = 0.1
    dt = TR 
    
    # extract signals and save them in preproc/roi_signals
    #-----------------------------------------------------
    run_4d = concat_niimgs(file_names, ensure_ndim=4)
    signals, _issues, _info = ucr.extract_signals(droi, roi_prefix, run_4d, 
                                                        mask=mask, minvox=1)   
    arr_sig, labels = ucr._dict_signals_to_arr(signals)

    # rm previous directory
    suf.rm_and_create(dsig)
    np.savez(fn_sig, arr_sig, labels) 

    # construct matrix of counfounds
    #-----------------------------------------------------

    #--- get CSF
    csf_signals, csf_issues, csf_info = \
                    ucr.extract_signals(dcsf, csf_filename, run_4d)
    #- check csf_signals ok ?
    assert len(csf_signals) == 1 # only one signal in dict
    csf_arr = csf_signals[csf_signals.keys()[0]]
    assert csf_arr.shape == (VOLNB,)
    csf_arr = csf_arr.reshape(VOLNB, 1)
    csf_labs = ['csf']
    # standardized csf_arr
    
    #--- get MVT
    mvt_arr = ucr.extract_mvt_run(mvtfile, idx_run, VOLNB)
    mvt_labs = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz']
    # standardized mvt_arr
    
    #--- get cosine functions;
    bf_arr, (order_lf, order_hf) = ucr._create_bandpass_bf(VOLNB, dt, 
                                                        (low_freq, high_freq))
    assert order_lf > 0 and order_hf > 0, "orders off {}, {}".format(
                                                        order_lf, order_hf)
    lf_labs = ["lf{:02d}".format(idx) for idx in range(order_lf)]
    hf_labs = ["hf{:02d}".format(idx) for idx in range(order_hf)]
    if verbose:
        print("order lf : {}, order hf : {}".format(order_lf, order_hf))
    
    
    #--- put it together  
    counfounds = np.hstack((csf_arr, mvt_arr, bf_arr))
    counfounds_labs = csf_labs + mvt_labs + lf_labs + hf_labs
    if verbose:
       print("csf.shape {}, mvt.shape {}, bf.shape {}".format(
                     csf_arr.shape, mvt_arr.shape, bf_arr.shape)

    # filter and compute correlation
    #-----------------------------------------------------

    flt_arr_sig = ucr._R_proj(counfounds, arr_sig)


    


    # save correlation




