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
import nibabel as nb
from nilearn._utils import concat_niimgs
import nibabel as nib

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
    argv[3]: session number, e.g. 01 or nothing - defaults to all sessions
    argv[4]: name of workflow to take data from - defaults to preproc
    """

    Ntrs = 195
    Conv_factor = (np.pi / 180) * 50

    ##1st arg: 01 (sub01); 2nd argument: 01 (ses01; optional, otherwise does all sessions)
    subject = 'sub' + str(argv[1])

    scrub_ver = argv[2]

    if len(sys.argv) > 3:
        session_info = 'sess' + str(argv[3])
    else:
        session_info = 'sess*'

    if len(sys.argv) > 4:
        workflow_name_func = argv[4]
    else:
        workflow_name_func = 'preproc'

    if scrub_ver == 'v1':
        Thr = [.08, .185, .335, .645]
    elif scrub_ver == 'v2':
        Thr = [.063, .153, .282, .55]
    elif scrub_ver == 'v3':
        Thr = [.15, .15, .15, .15]

    Flag_bound = [-1, 2]

    mot_conds = ['NONE', 'LOW', 'MED', 'HIGH']
    Nmotion = len(mot_conds)
    simpace_dir = '/home/despo/simpace'
    data_dir = pjoin(simpace_dir, 'rename_files')

    save_dir = pjoin(simpace_dir, 'analyses/fd-flag/' + subject)

    # get list of sess directories
    subject_dir = pjoin(data_dir, subject)

    # sessions = glob(subject_dir + '/' + session_info)

    Nsessions = 13
    Flag_sum = np.zeros((Nsessions, Nmotion))

    # for isess, sess_dir in enumerate(sessions):
    for isess in range(1, Nsessions + 1):

        if isess < 10:
            sess_str = '0' + str(isess)
        else:
            sess_str = str(isess)

        sess_str = 'sess' + sess_str
        sess_dir = pjoin(subject_dir, sess_str)

        sess = os.path.basename(sess_dir)
        workflow_dir = pjoin(sess_dir, workflow_name_func)

        print "working on session " + sess_dir

        # setup the directory with the motion params
        
        # load the order of motion conditions across runs
        params_fname = sess_dir + '/' + subject + '_' + sess + '_params.json'
        with open(params_fname) as fparam:
            params = json.load(fparam)

        mot_order = np.array(params['motion'])
        Nruns = len(mot_order)
        # pdb.set_trace()
        # mot_order = mot_order[[0,2,3]]
        # run_n = [1,3,4]

        # if not Nruns == 4:
        # print 'WARNING: ' + str(Nruns) + ' detected for ' subject + ', ' + sess

        # load the motion file and reshape to nruns x nvols
        motion_dir = pjoin(workflow_dir, 'motion_params')

        for irun in range(1,Nruns+1):
            motion_file = glob(pjoin(motion_dir, 'MP*run0' + str(irun) + '.txt'))[0]
            temp = np.loadtxt(motion_file)
            if irun==1:
                Nvols = temp.shape[0]
                Nmot_params = temp.shape[1]
                mot = np.zeros((Nvols, Nmot_params, Nruns))
            mot[:, :, irun-1] = temp

        # pipeline_dir = glob( pjoin(workflow_dir, pipeline) )
        # print pipeline_dir

        # mask_file = glob(pjoin(workflow_dir, 'sess_mask', 'sess_mask.nii'))
        # print mask_file
        #
        # mask = nib.load(mask_file[0])
        # mask_arr = np.asarray(mask.get_data()).astype(bool)

        FD = np.zeros((Ntrs, Nruns))
        Flag = np.zeros((Ntrs, Nruns))
        Flag_runs = np.zeros((Ntrs, Nruns))
        # DVARS_runs = np.zeros((Ntrs, Nruns))

        for irun, motion_cond in enumerate(mot_order):

            print motion_cond

            Mot_curr = mot[:, :, irun]
            # print Mot_curr.shape[0]
            # print Mot_curr.shape[1]

            AbsMotDiff = abs(np.diff(Mot_curr, axis=0))

            Trans = AbsMotDiff[:, 0:3]
            Rot = AbsMotDiff[:, 3:6]
            FD[1:Ntrs, irun] = np.sum(Trans, axis=1) + np.sum(Conv_factor * Rot, axis=1)

            # pdb.set_trace()
            for imot, cond_name in enumerate(mot_conds):
                if motion_cond == cond_name:
                    motion_idx = imot

            for iTR in np.arange(1, Ntrs + 1):

                if FD[iTR - 1, irun] > Thr[motion_idx]:

                    for iTR_rem in np.arange(iTR + Flag_bound[0], iTR + Flag_bound[1] + 1):
                        # print iTR_rem

                        if iTR_rem < Ntrs + 1:
                            Flag[iTR_rem - 1, irun] = 1

                        # pdb.set_trace()

            # Flag[:,irun] = FD[:,irun] > Thr[motion_idx]

            # FD_runs[:,motion_idx] = FD[:,irun]
            Flag_runs[:, motion_idx] = Flag[:, irun]

        # pdb.set_trace()
        Flag_sum[isess - 1, :] = np.sum(Flag_runs, axis=0)
        np.savetxt(pjoin(save_dir, sess_str + '_FD_flag_' + scrub_ver + '.txt'), Flag)

    print Flag_sum
    print Flag_sum / Ntrs
    print Ntrs - Flag_sum

    a = np.ma.array(Ntrs - Flag_sum, mask=False)
    a.mask[11, :] = True
    print np.mean(a, axis=0)


# pdb.set_trace()


if __name__ == '__main__':
    main()
