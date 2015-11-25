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
import dicom
import shutil
from six import string_types
import nibabel as nb


#-----------------------------------------------------------------------------
# Main Script
#-----------------------------------------------------------------------------
DESPO_SIMPACE_DIR = '/home/despo/simpace/rename_files/'
SUB_NUM = 1

#options
LW = 4

def main(argv = sys.argv):
    """
    argv[1]: subject number, e.g. 01
    argv[2]: session number, e.g. 01 or nothing - defaults to all sessions
    argv[3]: name of workflow to take data from - defaults to preproc
    """

    Conv_factor = (np.pi/180)*50
    print Conv_factor

    ##1st arg: 01 (sub01); 2nd argument: 01 (ses01; optional, otherwise does all sessions)
    subject = 'sub' + str(argv[1])

    if len(sys.argv) > 2:
        session_info = 'sess' + str(argv[2])
    else: 
        session_info = 'sess*'

    if len(sys.argv) > 3:
        workflow_name_func = argv[3]
    else:
        workflow_name_func = 'preproc'

    simpace_dir = '/home/despo/simpace'
    data_dir = pjoin(simpace_dir, 'rename_files')


    #get list of sess directories
    subject_dir = pjoin(data_dir, subject)

    sessions = glob(subject_dir + '/' + session_info)
    
    for isess, sess_dir in enumerate(sessions):

        sess = os.path.basename(sess_dir)
        workflow_dir = pjoin(sess_dir, workflow_name_func)

        print "working on session " + session
        
        #setup the directory with the motion params
        motion_dir = pjoin(workflow_dir, 'motion_params')
        motion_file = motion_dir + '/MP.txt'
        
        #load the order of motion conditions across runs
        params_fname = sess_dir + '/' + subject + '_' + sess + '_params' 
        with open(params_fname) as fparam:
            params = json.load(fparam)

        mot_order = np.array(params['motion'])
        Nruns = len(mot_order)
        if not Nruns == 4:
            print 'WARNING: ' str(Nruns) + ' detected for ' subject + ', ' + sess
        
        #load the motion file and reshape to nruns x nvols
        mot = np.loadtxt(motion_file)
        Nvols = mot.shape[0]/Nruns.
        Nmot_params = mot.shape[1]
        mot = np.reshape(mot, (Nruns,Nvols,Nmot_params))


        for irun, motion_cond in enumerate(mot_order):

            Mot_curr = mot[irun,:,:]
            print Mot_curr.shape[0]
            print Mot_curr.shape[1]

            AbsMotDiff = abs(np.diff(Mot_curr, axis = 0))

            Trans = AbsMotDiff[:,0:3]
            Rot = AbsMotD
            FD = np.sum(Trans, axis = 1) + np.sum(Conv_factor*Rot, axis = 1)
            print FD.shape[0]

           

    #     #set up plotting, for the xyz motion 
    #     fig = plt.figure()
    #     fig.set_size_inches(12,10,forward=True)
    #     mot_conds = ['NONE','LOW','MED','HIGH']
    #     for motidx, motcond in enumerate(mot_conds):
            
    #         #find the run in this sess that corresponds to motcond
    #         run = np.where(mot_order==motcond)[0][0]
            
    #         ax1 = fig.add_subplot(3,1,1)
    #         ax1.plot(mot[run,:,0], label=motcond,lw=LW)

    #         ax2 = fig.add_subplot(3,1,2)
    #         ax2.plot(mot[run,:,1], label=motcond,lw=LW)

    #         ax3 = fig.add_subplot(3,1,3)
    #         ax3.plot(mot[run,:,2], label=motcond,lw=LW)

    #     ax1.legend(loc='best')
    #     ax1.set_xticklabels([])
    #     ax2.set_xticklabels([])
    #     ax1.set_ylabel('X-axis translation (mm)', size=14)
    #     ax2.set_ylabel('Y-axis translation (mm)', size=14)
    #     ax3.set_ylabel('Z-axis translation (mm)', size=14)
    #     ax3.set_xlabel('Volume', size=16)
    #     fig.subplots_adjust(hspace = 0.15)
    #     fig.show()
    #     fig_fname = pjoin(sess, 'sub%02d_sess%02d_motion.pdf' %(int(sub),sessidx+1))
    #     fig.savefig(fig_fname)
    #     1/0

        

if __name__ == '__main__':
    main()
