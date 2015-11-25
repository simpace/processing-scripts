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
#import nipype.interfaces.freesurfer 

matlab.MatlabCommand.set_default_paths('/usr/local/matlab-tools/spm/spm8')

#1st arg: 01 (sub01); 2nd argument: 01 (ses01; optional, otherwise does all sessions)
subject = 'sub' + str(sys.argv[1])

if len(sys.argv) > 2:
	session_info = 'sess' + str(sys.argv[2])
else: 
	session_info = 'sess*'

#setup directories
#ASSUMES DIRECTORY STRUCTURE: data_dir/sub01/ses01/... and data_dir/sub01/anatomical/*.dcm
#ASSUMES atlas ROIs are in: simpace_dir/templates/atlas_name/roi_files (i.e. atlas_name='aal')
simpace_dir = '/home/despo/simpace'
data_dir = pjoin(simpace_dir, 'rename_files')
templates_dir = pjoin(simpace_dir, 'templates')

workflow_name_anatomical = 'anatomical_proc'
workflow_name_func = 'preproc'

nifti_dir = '/niftis/'
nii_ext = '/*.nii'
dicom_ext = '/*.dcm'
atlas_name = 'aal'
atlas_roi_files_dir = pjoin(templates_dir, atlas_name, 'roi_files')
atlas_roi_files_naming = atlas_name + '_*.nii'
		
#specific directories
subject_dir = pjoin(data_dir, subject)

anat_dir = pjoin(subject_dir, 'anatomical')
anat_proc_dir = pjoin(anat_dir, workflow_name_anatomical)
anat_file = glob(anat_dir + nifti_dir + nii_ext)

sessions = glob(subject_dir + '/' + session_info)
print subject_dir
print session_info
print sessions

#Run anatomical dcm --> nifti if anatomical_nifti dir doesn't exist
if not os.path.isdir(anat_proc_dir):
	TO_run_anatomical = True

	#template brain
	template_file = templates_dir + '/mni/avg152T1.nii'

	#aal files
	atlas_roi_files = glob(pjoin(atlas_roi_files_dir, atlas_roi_files_naming))
	Nrois = len(atlas_roi_files)
	print Nrois

else:
	TO_run_anatomical = False




if TO_run_anatomical:


	# anat_convert = pe.Node(interface=spmu.DicomImport(in_files=anat_files), name = 'anat_convert')

	segment = pe.Node(interface=spm.Segment(), name = 'segment')
	segment.inputs.data = anat_file
	segment.inputs.csf_output_type = [False,False,True]
	segment.inputs.gm_output_type = [False,False,True]
	segment.inputs.wm_output_type = [False,False,True]

	coreg_mni_to_native = pe.Node(interface=spm.Coregister(), name = 'coreg_mni_to_native')
	coreg_mni_to_native.inputs.target = anat_file[0]
	coreg_mni_to_native.inputs.source = template_file
	coreg_mni_to_native.inputs.apply_to_files = atlas_roi_files

	wm_file = pe.Node(interface=utility.Rename(), name = 'wm_file')
	wm_file.inputs.format_string = 'wm_anatres.nii'

	gm_file = pe.Node(interface=utility.Rename(), name = 'gm_file')
	gm_file.inputs.format_string = 'gm_anatres.nii'

	csf_file = pe.Node(interface=utility.Rename(), name = 'csf_file')
	csf_file.inputs.format_string = 'csf_anatres.nii'

	merge = pe.Node(interface=utility.Merge(2), name= 'merge')

	workflow_anatomical = pe.Workflow(name=workflow_name_anatomical)
	workflow_anatomical.base_dir = anat_dir

	workflow_anatomical.connect([ (segment, gm_file, [('native_gm_image', 'in_file')]),
		(segment, wm_file, [('native_wm_image', 'in_file')]),
		(segment, csf_file, [('native_csf_image', 'in_file')]),
		(coreg_mni_to_native, merge, [('coregistered_files','in1')]), 
		(coreg_mni_to_native, merge, [('coregistered_source', 'in2')]),
		])

	
	# coreg_dir = pjoin(anat_proc_dir, 'coreg_mni_to_native')
	# if not os.path.isdir(coreg_dir):
	# 	os.mkdir(anat_proc_dir)
	# 	os.mkdir(coreg_dir)

	# out_files = coreg_mni_to_native.run()
	# print out_files.coregistered_files


	# shutil.move(out_files.coregistered_source, coreg_dir + '/template_in_native.nii')
	# shutil.move(out_files.coregistered_files, coreg_dir )

	workflow_anatomical.run()

	##PUT IN MERGE OF SEGMENTED FILES + output of coreg_mni (template brain in native space)

	#merges takes input from segment --> rename above
	#	(segment, merge, [('native_gm_image','in1')]),
	#	(segment, merge, [('native_wm_image','in2')]),
	#	(segment, merge, [('native_csf_image','in3')]),
	#(coreg_mni, merge, [('coregistered_files','in4')]), #

if len(sessions) > 0:

	atlas_files = glob(anat_proc_dir + '/coreg_mni_to_native/r' + atlas_roi_files_naming)
	Nrois = len(atlas_files)

	files_to_coreg_to_native = [ glob(anat_proc_dir + '/gm_file' + nii_ext)[0],
		glob(anat_proc_dir + '/wm_file' + nii_ext)[0],
		glob(anat_proc_dir + '/csf_file' + nii_ext)[0] ]
	files_to_coreg_to_native.extend( atlas_files )
	print len(files_to_coreg_to_native)

		

for isess, sess_dir in enumerate(sessions):

	sess = os.path.basename(sess_dir)
	workflow_dir = pjoin(sess_dir, workflow_name_func)

	#grab dicom input files
	run_folders = subject + '_' + sess + '_run*'
	input_files = glob(sess_dir + '/' + run_folders + nifti_dir + nii_ext)

	
	#for fixing nipype file naming bug after smoothing 'sa' --> 'sra'
	out_dir = workflow_dir + '/smooth/'
	prefix_toremove = 'sa'
	prefix_toadd = 'sra'
	suffix = '*.nii'

	stc = pe.Node(interface=spm.SliceTiming(), name='stc')

	stc.inputs.in_files = input_files
	stc.inputs.num_slices = 32
	stc.inputs.time_repetition = 2.0
	stc.inputs.time_acquisition = 2. - 2./32
	stc.inputs.slice_order = range(32,0,-1) #prob not correct (clg-slices are descending)
	stc.inputs.ref_slice = 16 #close to middle slice in time

	realign = pe.Node(interface=spm.Realign(), name='realign')
	realign.inputs.register_to_mean = True

	smooth = pe.Node(interface=spm.Smooth(), name = 'smooth')
	smooth.inputs.fwhm = [6, 6, 6]

	coreg = pe.Node(interface=spm.Coregister(), name = 'coreg')
	coreg.inputs.source = anat_file
	coreg.inputs.apply_to_files = files_to_coreg_to_native #output of segment & roi_files from convert_mni_to_native

	select_gm = pe.Node(interface=utility.Select(), name = 'select_gm')
	select_gm.inputs.set(index=[0]) #selects gm

	select_wm = pe.Node(interface=utility.Select(), name = 'select_wm')
	select_wm.inputs.set(index=[1]) #selects wm

	select_csf = pe.Node(interface=utility.Select(), name = 'select_csf')
	select_csf.inputs.set(index=[2]) #selects csf

	# split & rename for aal outputs from coreg 4:(99+4)
	select_roi_files = pe.Node(interface=utility.Select(), name = 'select_roi_files')
	select_roi_files.inputs.set(index=range(3,3+Nrois))

	motion_params = pe.Node(interface=utility.Rename(), name = 'motion_params')
	motion_params.inputs.format_string = 'MP.txt'

	anat_in_func_res = pe.Node(interface=utility.Rename(), name = 'anat_in_func_res')
	anat_in_func_res.inputs.format_string = 'T1_func_res.nii'

	csf_in_func_res = pe.Node(interface=utility.Rename(), name = 'csf_in_func_res')
	csf_in_func_res.inputs.format_string = 'csf_func_res.nii'

	wm_in_func_res = pe.Node(interface=utility.Rename(), name = 'wm_in_func_res')
	wm_in_func_res.inputs.format_string = 'wm_func_res.nii'

	gm_in_func_res = pe.Node(interface=utility.Rename(), name = 'gm_in_func_res')
	gm_in_func_res.inputs.format_string = 'gm_func_res.nii'

	workflow_func = pe.Workflow(name=workflow_name_func)
	workflow_func.base_dir = sess_dir

	workflow_func.connect([ (stc, realign, [('timecorrected_files','in_files')]),
		(realign, smooth, [('realigned_files','in_files')]),
		(realign, coreg, [('mean_image','target')]),
		(coreg, select_csf, [('coregistered_files', 'inlist')]),
		(coreg, select_wm, [('coregistered_files', 'inlist')]),	
		(coreg, select_gm, [('coregistered_files', 'inlist')]),	
		(coreg, select_roi_files, [('coregistered_files', 'inlist')]),
		(coreg, anat_in_func_res, [('coregistered_source', 'in_file')]),
		(select_csf, csf_in_func_res, [('out', 'in_file')]), # NEED?
		(select_wm, wm_in_func_res, [('out', 'in_file')]),
		(select_gm, gm_in_func_res, [('out', 'in_file')]),
		(realign, motion_params, [('realignment_parameters', 'in_file')])
		])

	workflow_func.run()

	os.system('python Pmaps_v2.py ' + subject + ' ' + sess + ' ' + workflow_name_func)


	# FIX nipype file naming bug
	# full_file_name = glob(out_dir + prefix_toremove + suffix)
	# base_file_name = os.path.basename(full_file_name[0])
	# new_file_name = out_dir + prefix_toadd + base_file_name[len(prefix_toremove):]
	# os.rename(full_file_name[0], new_file_name)

	os.system('gzip ' + workflow_dir + '/stc/a*.nii')
	os.system('gzip ' + workflow_dir + '/realign/ra*.nii')

