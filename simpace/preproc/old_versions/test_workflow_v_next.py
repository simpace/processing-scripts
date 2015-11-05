"""Test workflow modified by CLG from test_workflow_at.py in 
/home/despo/arielle/simpace_testing/ 
"""

import os
import sys
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

sessN = sys.argv[1]
print sessN

simpace_dir = '/home/despo/simpace/'
subject = 'subject_1_data'
subject_othername = 'sub01'
session = 'sess' + sessN
nifti_dir = '/data'
data_dir = simpace_dir + subject + '/rename_files/' + subject_othername + '/' + session
func_dirs = pjoin(data_dir,'sub*sess*run*') #anat_dir = pjoin(data_dir,'*mprage*')
nii_ext = '/*.nii'
dicom_ext = '/*.dcm'
aal_name = 'aal_*.nii'
workflow_name = 'preproc'
past_dir = simpace_dir + subject + '/rename_files/sub01/sess01/' + workflow_name
save_dir = data_dir
print past_dir

#load the dicoms
input_files = glob(func_dirs + nifti_dir + nii_ext)

#template brain
template_file = simpace_dir + '/templates/mni/avg152T1.nii'

#aal files
aal_dir = simpace_dir + 'templates/aal/roi_files/'
roi_files = glob(aal_dir + aal_name)
Nrois = len(roi_files)
print Nrois

#anatomical files from previous run for coregistration into current functional space
gm_file = glob(past_dir + '/segmentmaps/c1*.nii')
wm_file = glob(past_dir + '/segmentmaps/c2*.nii')
csf_file = glob(past_dir + '/segmentmaps/c3*.nii')
roi_files = glob(past_dir + '/coreg_template/' + nii_ext)
apply_coreg_files_to = gm_file + wm_file + csf_file + roi_files
#print apply_coreg_files_to

#for fixing nipype file naming bug after smoothing 'sa' --> 'sra'
out_dir = save_dir + '/' + workflow_name + '/smooth/'
prefix_toremove = 'sa'
prefix_toadd = 'sra'
suffix = '*.nii'


#setup the nodes

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

coreg = pe.Node(interface=spm.Coregister(), name = 'registered_files')
coreg.inputs.source = glob(past_dir + '/anat/converted_dicom' + nii_ext)
coreg.inputs.apply_to_files = apply_coreg_files_to #output of merge
print coreg.inputs.source

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

workflow = pe.Workflow(name=workflow_name)
workflow.base_dir = save_dir

workflow.connect([ (stc, realign, [('timecorrected_files','in_files')]),
	(realign, smooth, [('realigned_files','in_files')]),
	(realign, coreg, [('mean_image','target')]),
	(coreg, split_csf, [('coregistered_files', 'inlist')]),
	(coreg, split_wm, [('coregistered_files', 'inlist')]),	
	(coreg, split_gm, [('coregistered_files', 'inlist')]),	
	(coreg, split_roi_files, [('coregistered_files', 'inlist')]),
	(coreg, rename_anat, [('coregistered_source', 'in_file')]),
	(split_csf, rename_csf, [('out', 'in_file')]),
	(split_wm, rename_wm, [('out', 'in_file')]),
	(realign, rename_mps, [('realignment_parameters', 'in_file')])
	])


#write and plot workflow graph
workflow.write_graph('test_workflow.dot')
#graph = misc.imread(pjoin(workflow.base_dir,'preproc','test_workflow.dot.png'))
#plt.imshow(graph)
#plt.show()

#run the workflow
workflow.run()


full_file_name = glob(out_dir + prefix_toremove + suffix)
base_file_name = os.path.basename(full_file_name[0])
new_file_name = out_dir + prefix_toadd + base_file_name[len(prefix_toremove):]
os.rename(full_file_name[0], new_file_name )

