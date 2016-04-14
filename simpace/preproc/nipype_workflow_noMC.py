"""Test workflow modified by CLG from test_workflow_at.py in 
/home/despo/arielle/simpace_testing/ 
"""

import os, sys, shutil
from glob import glob
from os.path import join as pjoin
from matplotlib import pyplot as plt
from scipy import misc
import nipype.interfaces.spm.utils as spmu
import nipype.interfaces.spm as spm
import nipype.interfaces.matlab as matlab
import nipype.interfaces.utility as utility
import nipype.pipeline.engine as pe
import nipype.interfaces.fsl.maths as fslmaths
from nipype.interfaces.freesurfer import DICOMConvert
from nipype.interfaces.afni import To3D
import json
import numpy as np
from nipype.interfaces import afni
import pdb

# import nipype.interfaces.freesurfer

matlab.MatlabCommand.set_default_paths('/usr/local/matlab-tools/spm/spm8')

# 1st arg: 01 (sub01); 2nd argument: 01 (ses01; optional, otherwise does all sessions)
subject = 'sub' + str(sys.argv[1])

if len(sys.argv) > 2:
    session_info = 'sess' + str(sys.argv[2])
else:
    session_info = 'sess*'

# setup directories
# ASSUMES DIRECTORY STRUCTURE: data_dir/sub01/ses01/... and data_dir/sub01/anatomical/*.dcm
# ASSUMES atlas ROIs are in: simpace_dir/templates/atlas_name/roi_files (i.e. atlas_name='aal')
simpace_dir = '/home/despo/simpace'
data_dir = pjoin(simpace_dir, 'rename_files')
templates_dir = pjoin(simpace_dir, 'templates')

workflow_name_anatomical = 'anatomical_proc'
workflow_name_func = 'preproc_noMC'

COREG_condition = 'NONE'
to_reg_atlas = False

nifti_dir = 'niftis'
nii_ext = '*.nii'
atlas_name = ['aal', 'power']

# atlas_roi_files_dir = pjoin(templates_dir, atlas_name, 'roi_files')
# atlas_roi_files_naming = atlas_name + '_*.nii'

# specific directories
subject_dir = pjoin(data_dir, subject)

anat_dir = pjoin(subject_dir, 'anatomical')
anat_proc_dir = pjoin(anat_dir, workflow_name_anatomical)
anat_file = glob(pjoin(anat_dir, nifti_dir, nii_ext))[0]

sessions = glob(subject_dir + '/' + session_info)
print subject_dir
print session_info
print sessions

if len(sessions) > 0:

    files_to_coreg_to_native = [glob(pjoin(anat_proc_dir, 'gm_file', nii_ext))[0],
                                glob(pjoin(anat_proc_dir, 'wm_file', nii_ext))[0],
                                glob(pjoin(anat_proc_dir, 'csf_file', nii_ext))[0]]

    if to_reg_atlas:
        coreg_dir = pjoin(anat_proc_dir, 'coreg_mni_to_native')

        # atlas files
        for iatlas, atlas in enumerate(atlas_name):
            atlas_roi_files_naming = atlas + nii_ext

            temp_list = glob(pjoin(coreg_dir, atlas_roi_files_naming))
            if iatlas == 0:
                atlas_files = temp_list
            else:
                atlas_files = atlas_files + temp_list

        Nrois = len(atlas_files)
        files_to_coreg_to_native.extend(atlas_files)

    print len(files_to_coreg_to_native)

for isess, sess_dir in enumerate(sessions):
    sess = os.path.basename(sess_dir)
    workflow_dir = pjoin(sess_dir, workflow_name_func)

    # # grab input files
    # run_folders = subject + '_' + sess + '_run*'
    # input_files = glob(sess_dir + '/' + run_folders + nifti_dir + nii_ext)
    # input_files.sort()

    run_str = subject + '_' + sess + '_run*'
    run_folders = glob(pjoin(sess_dir, run_str))
    input_files = []
    for irun in run_folders:
        input_files.append(glob(pjoin(irun, nifti_dir, os.path.basename(irun)) + nii_ext[1:len(nii_ext)])[0])
    input_files.sort()

    # load the order of motion conditions across runs
    params_fname = sess_dir + '/' + subject + '_' + sess + '_params.json'
    with open(params_fname) as fparam:
        params = json.load(fparam)

    mot_order = np.array(params['motion'])
    run_idx = np.where(mot_order == COREG_condition)[0][0]
    Ncond = len(mot_order)

    stc = pe.Node(interface=spm.SliceTiming(), name='stc')
    stc.inputs.in_files = input_files[run_idx]
    stc.inputs.num_slices = 32
    stc.inputs.time_repetition = 2.0
    stc.inputs.time_acquisition = 2. - 2. / 32
    stc.inputs.slice_order = range(32, 0, -1)  # prob not correct (clg-slices are descending)
    stc.inputs.ref_slice = 16  # close to middle slice in time

    smooth = pe.Node(interface=spm.Smooth(), name='smooth')
    smooth.inputs.fwhm = [6, 6, 6]

    select_smooth = pe.Node(interface=utility.Select(), name='select_smooth')
    select_smooth.inputs.set(index=[run_idx * (len(input_files) / Ncond)])  # selects image to coreg anatomical to

    coreg = pe.Node(interface=spm.Coregister(), name='coreg')
    coreg.inputs.source = anat_file
    coreg.inputs.apply_to_files = files_to_coreg_to_native  # output of segment & roi_files from convert_mni_to_native

    select_gm = pe.Node(interface=utility.Select(), name='select_gm')
    select_gm.inputs.set(index=[0])  # selects gm

    select_wm = pe.Node(interface=utility.Select(), name='select_wm')
    select_wm.inputs.set(index=[1])  # selects wm

    select_csf = pe.Node(interface=utility.Select(), name='select_csf')
    select_csf.inputs.set(index=[2])  # selects csf

    # split & rename for aal outputs from coreg 4:(99+4)
    # select_roi_files = pe.Node(interface=utility.Select(), name='select_roi_files')
    # select_roi_files.inputs.set(index=range(3, 3 + Nrois))

    anat_in_func_res = pe.Node(interface=utility.Rename(), name='anat_in_func_res')
    anat_in_func_res.inputs.format_string = 'T1_func_res.nii'

    csf_in_func_res = pe.Node(interface=utility.Rename(), name='csf_in_func_res')
    csf_in_func_res.inputs.format_string = 'csf_func_res.nii'

    wm_in_func_res = pe.Node(interface=utility.Rename(), name='wm_in_func_res')
    wm_in_func_res.inputs.format_string = 'wm_func_res.nii'

    gm_in_func_res = pe.Node(interface=utility.Rename(), name='gm_in_func_res')
    gm_in_func_res.inputs.format_string = 'gm_func_res.nii'

    workflow_func = pe.Workflow(name=workflow_name_func)
    workflow_func.base_dir = sess_dir

                           # (smooth, select_smooth, [('smoothed_files', 'inlist')]),

    workflow_func.connect([(stc, smooth, [('timecorrected_files', 'in_files')])])

    workflow_func.run()

    os.system('gzip ' + workflow_dir + '/stc/a*.nii')
