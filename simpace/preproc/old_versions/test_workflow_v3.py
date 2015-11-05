"""Test workflow modified by CLG from test_workflow_at.py in 
/home/despo/arielle/simpace_testing/ 
"""

import os
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

#setup directories
simpace_dir = '/home/despo/simpace/'
subject = 'subject_1_data'
subject_othername = 'sub01'
session = 'sess01'
nifti_dir = '/data'
data_dir = simpace_dir + subject + '/rename_files/' + subject_othername + '/' + session
anat_dir = simpace_dir + subject + '/anatomical_may_12_2015/*MPRAGE*' #PREVIOUSLY: '/temp_anatomical/*mprage*'
func_dirs = pjoin(data_dir,'sub*sess*run*') #
nii_ext = '/*.nii'
dicom_ext = '/*.dcm'
aal_name = 'aal_*.nii'
workflow_name = 'preproc'
save_dir = data_dir
print func_dirs

#load the dicoms
input_files = glob(func_dirs + nifti_dir + nii_ext)
anat_files = glob(anat_dir + dicom_ext)

#template brain
template_file = simpace_dir + '/templates/mni/avg152T1.nii'

#aal files
aal_dir = simpace_dir + 'templates/aal/roi_files/'
roi_files = glob(aal_dir + aal_name)
Nrois = len(roi_files)
print Nrois

#for fixing nipype file naming bug after smoothing 'sa' --> 'sra'
out_dir = save_dir + '/' + workflow_name + '/smooth/'
prefix_toremove = 'sa'
prefix_toadd = 'sra'
suffix = '*.nii'

#setup the nodes

anat_convert = pe.Node(interface=spmu.DicomImport(in_files=anat_files), name='anat')

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

segment = pe.Node(interface=spm.Segment(), name = 'segmentmaps')
segment.inputs.csf_output_type = [False,False,True]
segment.inputs.gm_output_type = [False,False,True]
segment.inputs.wm_output_type = [False,False,True]

merge = pe.Node(interface=utility.Merge(4), name= 'merge_files_for_epicoreg')

coreg = pe.Node(interface=spm.Coregister(), name = 'registered_files')

coreg_mni = pe.Node(interface=spm.Coregister(source = template_file), name = 'coreg_template')
coreg_mni.inputs.source = template_file
coreg_mni.inputs.apply_to_files = roi_files

split_csf = pe.Node(interface=utility.Select(), name = 'split_csfmap')
split_csf.inputs.set(index=[2]) #selects csf

split_wm = pe.Node(interface=utility.Select(), name = 'split_wmmap')
split_wm.inputs.set(index=[1]) #selects wm

split_gm = pe.Node(interface=utility.Select(), name = 'split_gmmap')
split_gm.inputs.set(index=[0]) #selects gm

#ADD split & rename for aal outputs from coreg 4:(99+4)
split_roi_files = pe.Node(interface=utility.Select(), name = 'split_roi_files')
split_roi_files.inputs.set(index=range(3,3+Nrois))

rename_wm = pe.Node(interface=utility.Rename(), name = 'rename_wm')
rename_wm.inputs.format_string = 'wm_map.nii'

rename_csf = pe.Node(interface=utility.Rename(), name = 'rename_csf')
rename_csf.inputs.format_string = 'csf_map.nii'

rename_mps = pe.Node(interface=utility.Rename(), name = 'rename_mps')
rename_mps.inputs.format_string = 'MP.txt'

rename_anat = pe.Node(interface=utility.Rename(), name = 'rename_anat_coreg')
rename_anat.inputs.format_string = 'T1_coreg.nii'

workflow = pe.Workflow(name='preproc')
workflow.base_dir = save_dir

workflow.connect([ (stc, realign, [('timecorrected_files','in_files')]),
	(realign, smooth, [('realigned_files','in_files')]),
	(anat_convert, segment, [('out_files','data')]),
	(anat_convert, coreg_mni, [('out_files','target')]),
	(anat_convert, coreg, [('out_files','source')]),
	(realign, coreg, [('mean_image','target')]),
	(segment, merge, [('native_gm_image','in1')]),
	(segment, merge, [('native_wm_image','in2')]),
	(segment, merge, [('native_csf_image','in3')]),
	(coreg_mni, merge, [('coregistered_files','in4')]), #
	(merge, coreg, [('out', 'apply_to_files')]), #coreg anatomical subject to functional space
	(coreg, split_csf, [('coregistered_files', 'inlist')]),
	(coreg, split_wm, [('coregistered_files', 'inlist')]),	
	(coreg, split_gm, [('coregistered_files', 'inlist')]),	
	(coreg, split_roi_files, [('coregistered_files', 'inlist')]),
	(coreg, rename_anat, [('coregistered_source', 'in_file')]),
	(split_csf, rename_csf, [('out', 'in_file')]),
	(split_wm, rename_wm, [('out', 'in_file')]),
	(realign, rename_mps, [('realignment_parameters', 'in_file')])
	])


	# (split_csf, threshold_csf_99, [('out', 'in_file')]),
	# (split_wm, threshold_wm_99, [('out', 'in_file')]),
#	(coreg, threshold, [('coregistered_files','in_file')])


#write and plot workflow graph
workflow.write_graph('test_workflow.dot')
#graph = misc.imread(pjoin(workflow.base_dir,'preproc','test_workflow.dot.png'))
#plt.imshow(graph)
#plt.show()

#run the workflow
workflow.run()
#  parameter plugin : sge


# FIX nipype file naming bug
full_file_name = glob(out_dir + prefix_toremove + suffix)
base_file_name = os.path.basename(full_file_name[0])
new_file_name = out_dir + prefix_toadd + base_file_name[len(prefix_toremove):]
os.rename(full_file_name[0], new_file_name )