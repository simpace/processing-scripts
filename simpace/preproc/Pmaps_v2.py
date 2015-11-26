import os, glob, sys
from os.path import join as pjoin
import nipype.interfaces.fsl.maths as fslmaths
import nipype.interfaces.fsl as fsl
import numpy as np
import SimpleITK as sitk
import nibabel as nb

# subject = 'sub' + str(sys.argv[1])
# session = 'sess' + str(sys.argv[2])
subject = sys.argv[1]
session = sys.argv[2]
if len(sys.argv) > 3:
	workflow_name = sys.argv[3]
else:
	workflow_name = 'preproc'

simpace_dir = '/home/despo/simpace'
data_dir = pjoin(simpace_dir, 'rename_files')

sess_dir = pjoin(data_dir, subject, session, workflow_name)
file_ext = '.nii'
gzip_ext = '.nii.gz'

print 'eroding masks for session: ' + sess_dir

folder_base = '_in_func_res'
masks = ['wm', 'csf']

Nvoxels_min = 20 #minimum # of voxels required after mask erosion

probability_thrs = np.arange(.80, .99, .01)

imagemaths = fsl.ImageMaths()

erodeITK = sitk.BinaryErodeImageFilter()
erodeITK.SetKernelRadius(1)
erodeITK.SetKernelType(3)
erodeITK.SetForegroundValue(1)

for imask in masks:

	mask_path = pjoin(sess_dir, imask + folder_base)
	file_path = mask_path + '/' + imask + '*' + file_ext
	mask_file = glob.glob(file_path)[0]
	mask_base = os.path.splitext(os.path.basename(mask_file))[0]

	final_mask_file = mask_path + '/' + mask_base + '_final' + gzip_ext

	for ithresh in probability_thrs:

		new_file_base = mask_path + '/' + mask_base + '_' + str(ithresh)
		thresh_file = new_file_base + gzip_ext
		erode_file = new_file_base + '_erode' + gzip_ext
		
		imagemaths.inputs.in_file = mask_file
		imagemaths.inputs.op_string = '-thr ' + str(ithresh) + ' -bin' #Threshold & binarize
		imagemaths.inputs.out_file = thresh_file
		imagemaths.run()

		img = sitk.ReadImage(thresh_file)
		img = sitk.Cast(img, sitk.sitkInt8)

		out = erodeITK.Execute(img)

		if imask == 'wm': #erode twice
			out_again = erodeITK.Execute(out)
			sitk.WriteImage(out_again,erode_file)

		elif imask == 'csf': #erode once
			sitk.WriteImage(out,erode_file)
	
		output_mask = nb.load(erode_file).get_data()
		total = np.sum(output_mask)
	
		if total >= Nvoxels_min:
			print str(total) + ' voxels in ' + erode_file

			if os.path.islink(final_mask_file):
				os.remove(final_mask_file)

			os.symlink( os.path.splitext(os.path.basename(erode_file))[0] + '.gz', final_mask_file)
