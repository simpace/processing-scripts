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
from nipype.algorithms import metrics

import pdb, traceback


# -----------------------------------------------------------------------------
# Main Script
# -----------------------------------------------------------------------------


def main(argv=sys.argv):
    """
    argv[1]: subject number, e.g. 01
    argv[2]: name of subfolder in workflow to take data from (PIPELINE) - no default, e.g. NR_4
    argv[3]: scrubbing info, False=No scrubbing, other value indicates scrubbing version, e.g. v2
    argv[4]: atlas, defaults to aal
    argv[5]: session number, e.g. 01 or nothing - defaults to all sessions
    argv[6]: name of workflow to take data from - defaults to preproc
    """

    subject = 'sub' + str(argv[1])

    pipeline = argv[2]

    scrub_info = argv[3]
    if 'False' in scrub_info:
        TOscrub = False
    else:
        TOscrub = True
        scrub_ver = scrub_info
    print TOscrub

    if len(sys.argv) > 4:
        roi_atlas = argv[4]
    else:
        roi_atlas = 'aal'

    if len(sys.argv) > 5:
        to_LFF = argv[5]
        if to_LFF in ['false', 'False']:
            to_LFF = False
        elif to_LFF in ['true', 'True']:
            to_LFF = True
    else:
        to_LFF = False

    if len(sys.argv) > 6:
        session_info = 'sess' + str(argv[6])
    else:
        session_info = 'sess*'

    if len(sys.argv) > 7:
        workflow_name_func = argv[7]
    else:
        workflow_name_func = 'preproc'

    # mot_conds = ['NONE', 'LOW', 'MED', 'HIGH']
    simpace_dir = '/home/despo/simpace'
    data_dir = pjoin(simpace_dir, 'rename_files')

    save_dir = pjoin(simpace_dir, 'analyses/roi-pair-distance/' + subject)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    # get list of sess directories
    subject_dir = pjoin(data_dir, subject)
    sessions = glob(subject_dir + '/' + session_info)

    Nsessions = len(sessions)
    Nruns = 4
    roi_threshold = .5
    Ntrs = 195

    signal_file_append = '*' + roi_atlas #+ '.npz'
    roi_sub_dir = 'coreg'

    if to_LFF:
        # signal_file_append += '_LFF'
        save_name_append_base = '_LFF'
        data_type = 'arr_sig_f'
    else:
        save_name_append_base = ''
        data_type = 'arr_sig'

    if TOscrub:
        save_name_append_base += '_SCRUB_' + scrub_ver + '_'

    Nrois_sess = np.zeros((Nsessions, 1))

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

        signals_dir = pjoin(workflow_dir, 'extracted_signals')
        roi_dir = pjoin(workflow_dir, roi_sub_dir)

        roi_save_file = pjoin(save_dir, pipeline + '_' + sess_str + '_' + roi_atlas + '_ROI_')
        print roi_save_file

        params_fname = sess_dir + '/' + subject + '_' + sess + '_params.json'
        with open(params_fname) as fparam:
            params = json.load(fparam)

        mot_order = np.array(params['motion'])

        mask = nib.load(pjoin(workflow_dir, 'sess_mask/sess_mask.nii'))
        mask_arr = np.asarray(mask.get_data()).astype(bool)

        if TOscrub:
            sc_dir = '/home/despo/simpace/analyses/fd-flag/' + subject
            vol_idx = np.loadtxt(glob(pjoin(sc_dir, sess_str + '*' + scrub_ver + '*.txt'))[0])
        else:
            vol_idx = np.zeros((Ntrs, Nruns))

        print 'using # of volumes: ' + str(np.sum(vol_idx == 0, axis=0))

        # get corr matrices
        for irun in range(1, Nruns + 1):

            npz_file = glob(pjoin(signals_dir, pipeline + '*run0' + str(irun) + signal_file_append + '_LFF.npz'))[0]
            data = np.load(npz_file)

            roi_sig = np.transpose(data[data_type])  # roi_sig is Nrois x Ntrs
            orig_shape = roi_sig.shape

            tr_idx = vol_idx[:, irun - 1] == 0
            roi_sig = roi_sig[:, tr_idx]

            if not TOscrub:
                assert orig_shape == roi_sig.shape
                save_name_append = save_name_append_base
            else:
                save_name_append = save_name_append_base + str(np.sum(tr_idx, axis=0))

            corr_mtx = np.corrcoef(roi_sig)
            save_file_nm = pjoin(save_dir, pipeline + '_' + sess_str + '_run0' + str(irun) + '_' + \
                                 mot_order[irun - 1] + '_' + roi_atlas)
            corr_save_file = save_file_nm + '_CorrMtx' + save_name_append + '.txt'

            # print corr_save_file
            # pdb.set_trace()
            np.savetxt(corr_save_file, corr_mtx)

            if irun == Nruns:
                roi_names = data['labels_sig']
                Nrois = len(roi_names)

        # get ROI information for each session
        distance_file = roi_save_file + 'Distance.txt'
        if not os.path.exists(distance_file):

            Nrois_sess[isess - 1] = Nrois
            Coord = np.zeros((Nrois, 3))
            Dist = np.zeros((Nrois, Nrois))

            # get all ROI coordinates
            for iroi, roi_curr in enumerate(roi_names):
                roi_file = glob(pjoin(roi_dir, roi_curr + '*'))
                roi = nib.load(roi_file[0])

                # roi_arr = np.nonzero(np.asarray(roi.get_data().astype(float) > roi_threshold) == 1)
                roi_arr = np.asarray(roi.get_data().astype(float) > roi_threshold) == 1
                roi_mask = np.nonzero(np.logical_and(roi_arr, mask_arr) == 1)

                affine = roi.affine

                Coord[iroi, 0:3] = np.dot(affine, np.hstack((np.mean(roi_mask, axis=1), 1)))[:3]

            # if iroi==44 and isess==5:
            # 	pdb.set_trace()

            # compute ROI pair distances
            for iroi in range(0, Nrois):

                for iroipair in range(iroi + 1, Nrois):
                    Dist[iroi, iroipair] = np.sqrt(np.sum(np.square(Coord[iroi, :] - Coord[iroipair, :])))

            np.savetxt(distance_file, Dist)

        # pdb.set_trace()
        # np.savetxt( roi_save_file + 'Names_' + str(Nrois) + '.txt', roi_names.as_list )

        roi_name_file = roi_save_file + 'Names_' + str(Nrois) + '.txt'

        if not os.path.exists(roi_name_file):
            np.savetxt(roi_name_file, roi_names, delimiter="", fmt="%s")


if __name__ == '__main__':
    main()
