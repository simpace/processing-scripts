import os, glob, sys
from os.path import join as pjoin
import numpy as np
from glob import glob
import nibabel as nb


#1st arg: 01 (sub01); 2nd argument: 01 (ses01; optional, otherwise does all sessions)
subject = 'sub' + str(sys.argv[1])

if len(sys.argv) > 2:
	session_info = 'sess' + str(sys.argv[2])
else: 
	session_info = 'sess*'

if len(sys.argv) > 3:
    workflow_name_func = argv[3]
else:
    workflow_name_func = 'preproc'

Nruns = 4

func_dir = 'smooth'
func_file_ext = '*.nii'

#setup directories
#ASSUMES DIRECTORY STRUCTURE: data_dir/sub01/ses01/... and data_dir/sub01/anatomical/*.dcm
#ASSUMES atlas ROIs are in: simpace_dir/templates/atlas_name/roi_files (i.e. atlas_name='aal')
simpace_dir = '/home/despo/simpace'
data_dir = pjoin(simpace_dir, 'rename_files')


#specific directories
subject_dir = pjoin(data_dir, subject)
sessions = glob(subject_dir + '/' + session_info)



for isess, sess_dir in enumerate(sessions):

    sess = os.path.basename(sess_dir)
    workflow_dir = pjoin(sess_dir, workflow_name_func)

    print "working on session " + workflow_dir
    
    #IF regressing out motion
    #setup the directory with the motion params
    motion_dir = pjoin(workflow_dir, 'motion_params')
    motion_file = motion_dir + '/MP.txt'

	#load the motion file and reshape to nruns x nvols
    mot = np.loadtxt(motion_file)
    Nvols = mot.shape[0]/Nruns
    Nmot_params = mot.shape[1]
    mot = np.reshape(mot, (Nruns,Nvols,Nmot_params))

    #IF WM
    wm_file = glob( pjoin(workflow_dir, 'wm_in_func_res') + '/wm_func_res_final.nii.gz')

    #IF CSF
    csf_file = glob( pjoin(workflow_dir, 'csf_in_func_res') + '/csf_func_res_final.nii.gz')

    # print wm_file
    # print csf_file

    csf_data = nb.load(csf_file[0]).get_data()
    
    Nvoxels = csf_data.size
    print Nvoxels
    csf_data = csf_data.reshape(Nvoxels)

    wm_data = nb.load(wm_file[0]).get_data()
    wm_data = wm_data.reshape(Nvoxels)
	
    print np.amax(wm_data)
    print np.amin(wm_data)
    

    func_file_dir = pjoin(workflow_dir, func_dir)
    
    for irun in range(1, 2): #Nruns+1):

    	files = glob( func_file_dir + '/*run0' + str(irun) + func_file_ext )

    	for ifile, file_name in enumerate(files):

    		# print ifile
    		# print file_name
	    		
    		data_perfile = nb.load(file_name).get_data()

    		if ifile == 0:
    			data = np.zeros( (Nvoxels, Nvols) )
    			print data.shape[0]
    			print data.shape[1]

	    	data[:, ifile] = data_perfile.reshape(Nvoxels)



		#CUE regressors for run

        Nuis = np.empty((0,Nvols), int)
        print Nuis.shape

		#IF CSF
        csf = data[csf_data==1, :]
        #IF mean
        csf = np.mean(csf, axis=0)

        Nuis = np.vstack((Nuis, csf))

        #IF WM
        wm = data[wm_data==1, :]
        #IF mean
        wm = np.mean(wm, axis=0)

        Nuis = np.vstack((Nuis, wm))
        print Nuis.shape

        #IF MOTION
        Nuis = np.vstack((Nuis, np.transpose(mot[irun,:,:])))


        print Nuis.shape





	    	# print irun
	    	# print files[0]
	    	# print data.shape