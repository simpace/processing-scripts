"""Script to plot motion parameters across runs/conditions for each session
"""

import os, sys, stat
import numpy as np
import re
import json
from glob import glob
from os.path import join as pjoin
from matplotlib import pyplot as plt
import os.path as osp
import shutil
from six import string_types
from nilearn._utils import concat_niimgs
import nibabel as nib
from simpace.utils import _utils as ucr
# from utils import _utils as ucr


import pdb

# -----------------------------------------------------------------------------
# Main Script
# -----------------------------------------------------------------------------
DESPO_SIMPACE_DIR = '/home/despo/simpace/rename_files/'

# options
LW = 4


def main(argv=sys.argv):
    """
    argv[1]: subject number, e.g. 01
    argv[2]: name of subfolder in workflow to take data from (PIPELINE) - no default
    argv[3]: name of workflow to take data from - defaults to preproc
    argv[4]: to mask by gray matter, True or False - default is False
    argv[5]: to perform LFF - default is False (no filtering)
    argv[6]: session number, e.g. 01 or nothing - defaults to all sessions

    """

    Ntrs = 195
    Conv_factor = (np.pi / 180) * 50
    TR = 2

    subject = 'sub' + str(argv[1])
    pipeline = argv[2]
    pipeline_savenm = pipeline

    if len(sys.argv) > 3:
        workflow_name_func = argv[3]
    else:
        workflow_name_func = 'preproc'

    workflow_name_save = workflow_name_func + '_'

    if len(sys.argv) > 4:
        to_mask_gm = argv[4]
        if to_mask_gm in ['false', 'False']:
            to_mask_gm = False
        elif to_mask_gm in ['true', 'True']:
            to_mask_gm = True
    else:
        to_mask_gm = False

    if len(sys.argv) > 5:
        to_LFF = argv[5]
        if to_LFF in ['false', 'False']:
            to_LFF = False
        elif to_LFF in ['true', 'True']:
            to_LFF = True
    else:
        to_LFF = False

    if to_mask_gm:
        gm_thr = .5
        type_savenm = '_GM'
    else:
        type_savenm = ''

    if len(sys.argv) > 6:
        session_info = 'sess' + str(argv[6])
    else:
        session_info = 'sess*'

    if to_LFF:
        type_savenm += '_LFF'

    mot_conds = ['NONE', 'LOW', 'MED', 'HIGH']
    Nmotion = len(mot_conds)

    simpace_dir = '/home/despo/simpace'
    data_dir = pjoin(simpace_dir, 'rename_files')
    save_dir = pjoin(simpace_dir, 'analyses/fd-dvars/' + subject)

    # get list of session directories
    subject_dir = pjoin(data_dir, subject)
    sessions = glob(subject_dir + '/' + session_info)

    Nsessions = len(sessions)
    Nmot_params = 6

    correl_sess = np.zeros((Nsessions, Nmotion))
    dp_sess = np.zeros((Nsessions, Nmotion))

    # for isess, sess_dir in enumerate(sessions):
    for isess in range(1, Nsessions + 1):

        if isess < 10:
            sess_str = '0' + str(isess)
        else:
            sess_str = str(isess)

        sess_dir = pjoin(subject_dir, 'sess' + sess_str)

        sess = os.path.basename(sess_dir)
        workflow_dir = pjoin(sess_dir, workflow_name_func)

        print "working on session " + sess_dir

        # setup the directory with the motion params
        motion_dir = pjoin(workflow_dir, 'motion_params')
        motion_files = glob(pjoin(motion_dir, 'MP*.txt'))

        # load the order of motion conditions across runs
        params_fname = pjoin(sess_dir, subject + '_' + sess + '_params.json')
        with open(params_fname) as fparam:
            params = json.load(fparam)

        mot_order = np.array(params['motion'])
        Nruns = len(mot_order)

        # decide whether or not to compute FD - check if file exists
        fd_file = pjoin(save_dir, sess_str + '_' + workflow_name_save + pipeline_savenm + '_FD.txt')
        if not os.path.exists(fd_file):
            to_compute_FD = True
        else:
            to_compute_FD = False

        if not os.path.exists(motion_files[0]):
            to_compute_FD = False

        if to_compute_FD:
            motion = np.zeros((Ntrs, Nmot_params, Nruns))
            for ifile, motion_file in enumerate(motion_files):
                motion[:, :, ifile] = np.loadtxt(motion_file)
            Nvols = motion.shape[0]
            if not Nvols == Ntrs:
                Nruns = motion.shape[0] / (Ntrs)
                print '*****only ' + str(Nruns) + ' RUN detected'

        pipeline_dir = glob(pjoin(workflow_dir, pipeline))[0]
        print pipeline_dir

        # load in session mask
        mask_file = glob(pjoin(workflow_dir, 'sess_mask', 'sess_mask.nii'))
        print mask_file
        mask = nib.load(mask_file[0])
        mask_arr = np.asarray(mask.get_data()).astype(bool)

        # additionally mask gray matter
        if to_mask_gm:
            gm_file = glob(pjoin(workflow_dir, 'gm_in_func_res', 'gm_func_res.nii'))
            gm = nib.load(gm_file[0])
            gm_mask = np.asarray(gm.get_data() > gm_thr)
            mask_total = np.logical_and(mask_arr, gm_mask)
        else:
            mask_total = mask_arr

        # initialize variables for FD-DVARS relationships
        correl = np.zeros((Nruns, 1))
        dp = np.zeros((Nruns, 1))

        DVARS = np.zeros((Ntrs-1, Nruns))
        DVARS_runs = np.zeros((Ntrs-1, Nruns))
        FD = np.zeros((Ntrs-1, Nruns))
        if to_compute_FD:
            FD_runs = np.zeros((Ntrs-1, Nruns))
        else:
            FD_runs = np.loadtxt(fd_file)

        for irun, motion_cond in enumerate(mot_order):
            print motion_cond

            if to_compute_FD:
                Mot_curr = motion[:, :, irun]

                AbsMotDiff = abs(np.diff(Mot_curr, axis=0))

                Trans = AbsMotDiff[:, 0:3]
                Rot = AbsMotDiff[:, 3:6]
                FD[:, irun] = np.sum(Trans, axis=1) + np.sum(Conv_factor * Rot, axis=1)

            else:
                junk = np.asarray(mot_conds, dtype='<U4') == motion_cond
                run_idx = np.where(junk)[0][0]
                FD[:, irun] = FD_runs[:, run_idx]

            run_str = '*run0' + str(irun + 1) + '*.nii.gz'  # indexing starts at 0

            print run_str
            if 'localWMreg' in pipeline:
                run_str = '*preproc_localWMreg*' + run_str + '*'

            file = glob(pjoin(pipeline_dir, run_str))
            if len(file) > 1:
                print len(file) + 'files detected for run ' + str(irun+1)

            # load in data
            run_4d = nib.load(file[0])
            img_arr = np.asarray(run_4d.get_data())
            data = img_arr[mask_total]

            if to_LFF:
                bf_arr, labels = ucr.extract_bf(float("NaN"), 0.08, Ntrs, TR)
                data_re = np.transpose(data)
                arr_sig_f = ucr.R_proj(bf_arr, data_re)
                data = np.transpose(arr_sig_f)

            DVARS[:, irun] = np.sqrt(np.mean(np.square(np.diff(data, 1, axis=1)), axis=0))

            temp = np.corrcoef(FD[:, irun], DVARS[:, irun])
            correl[irun] = temp[0, 1]
            dp[irun] = np.dot(FD[:, irun], DVARS[:, irun])

        # sort measures by motion level: None - High
        for icond, motion_cond in enumerate(mot_conds):
            motion_cond = mot_conds[icond]
            print icond
            print motion_cond
            run_idx = np.where(mot_order == motion_cond)[0][0]

            correl_sess[isess - 1, icond] = correl[run_idx]
            dp_sess[isess - 1, icond] = dp[run_idx]

            FD_runs[:, icond] = FD[:, run_idx]
            DVARS_runs[:, icond] = DVARS[:, run_idx]
            print run_idx

        print correl_sess
        print dp_sess

        if to_compute_FD:
            if not os.path.isdir(save_dir):
                os.mkdir(save_dir)
            np.savetxt(fd_file, FD_runs)

        dvars_file = pjoin(save_dir, sess_str + '_' + workflow_name_save + pipeline_savenm + '_DVARS' + type_savenm + '.txt')
        np.savetxt(dvars_file, DVARS_runs)

    out_file = pjoin(save_dir, workflow_name_save + pipeline_savenm + type_savenm)
    np.savetxt(out_file + '_correlation.txt', correl_sess)
    np.savetxt(out_file + '_dotproduct.txt', dp_sess)


if __name__ == '__main__':
    main()
