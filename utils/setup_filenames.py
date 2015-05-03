"""Script to setup dicom directories to standard format
"""

import os, sys, stat
import numpy as np
import re
import json
from glob import glob
from os.path import join as pjoin
import dicom
import shutil
import tempfile
import nipype.interfaces.spm.utils as spmu
import nipype.interfaces.matlab as matlab

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def isgroup_readable(filepath):
    st = os.stat(filepath)
    return bool(st.st_mode & stat.S_IRGRP)

def isgroup_writeable(filepath):
    st = os.stat(filepath)
    return bool(st.st_mode & stat.S_IWGRP)

def isgroup_executeable(filepath):
    st = os.stat(filepath)
    return bool(st.st_mode & stat.S_IXGRP)

def isuser_readable(filepath):
    st = os.stat(filepath)
    return bool(st.st_mode & stat.S_IRUSR)

def isuser_writeable(filepath):
    st = os.stat(filepath)
    return bool(st.st_mode & stat.S_IWUSR)

def isuser_executeable(filepath):
    st = os.stat(filepath)
    return bool(st.st_mode & stat.S_IXUSR)

def check_dir(filepath, checklist=['exists']):
    """
    check that directory has the attribute in checklist
    With no argument, check_dir verifies existence and
    is equivalent to os.path.isdir(filepath)
    parameters:
    -----------
    filename: string
        the filepath
    checklist: list of string
        with all the check to be performed
        ['exists', 'iswritable']
    returns:
    -------
    bool 
        a boolean if all checklist conditions are true
    """
    checkbool = True
    for chck in checklist:
        if chck == 'exists':
            checkbool = checkbool and os.path.isdir(filepath)
        elif chck == 'iswriteable':
            # try to write in that directory
            try:
                temp = tempfile.TemporaryFile(mode='w+b', dir=filepath)
                temp.close()
            except:
                checkbool = False

    return checkbool

def sort_ctime(file_list):
    """This function takes in a list of files and returns a sorted list of them 
    by their creation time"""
    
    files_sorted = sorted(file_list, key=os.path.getctime)

    return files_sorted


def get_nonorm(dcm_dirs):
    """This function reads in a list of directories with dicom files,
    loads a test dicom file from each and returns which is NOT prescan normalized
    by reading the dicom header"""

    #initialize a list for non prescan normed files
    no_psnorm = []

    #loop through the dicoms
    for dcmidx, dcm in enumerate (dcm_dirs):

        #load the first dicom in each directory)
        test_dcm = sorted(glob(pjoin(dcm,'IM-*.dcm')))[0]
        hdr = dicom.read_file(test_dcm)

        #test the ImageType field for 'NORM', if yes it's prescan normalized
        imtype = hdr.ImageType
        if 'NORM' not in imtype:
            no_psnorm.append(dcm)

    #TEST that there is only one file being returned
    if not len(no_psnorm) == 1:
        raise ValueError("more than one file matches the norm condition!")

    return no_psnorm


def sort_acqtime(dcm_dirs):
    """This function takes in a list of directories with dicom files,
    loads a test dicom file from each and returns them in order of acquisition"""

    #initialize an array to store the acquisition times for each dir
    acq_times = np.ones((len(dcm_dirs)))*np.nan

    #loop through the dicoms
    for dcmidx, dcm in enumerate (dcm_dirs):

        #load the first dicom in each directory)
        test_dcm = sorted(glob(pjoin(dcm,'IM-*.dcm')))[0]
        hdr = dicom.read_file(test_dcm)

        #get the acq time and store it in the array
        acq_times[dcmidx] = hdr.AcquisitionTime

    #sort the array on acq times, apply it to the dir names
    acq_order = acq_times.argsort()
    dcm_acqsort = np.array(dcm_dirs)[acq_order]

    return dcm_acqsort.tolist()
    

def parse_string(in_string, search_string):
    """Finds the index of everything up to a certain string element, returns all 
    portions of the string after this"""

    search = """.*%s(.*)""" %(search_string)
    match = re.search(search, in_string) #search the input for a search string
    new_string = match.groups() #find the elements after the search string
    #new_string = in_string[len(pre_search):]
    
    return new_string


def get_motorder(dcm_dirs):
    """This function takes in a list of ordered dicom directories and returns just the 
    motion condition relevant names"""

    #initialize a list of the order names
    mot_order = []

    #loop through the directory names
    for dcmidx, dcm in enumerate(dcm_dirs):
        
        #use the string parser to remove the initial dicom elements
        new_string = parse_string(dcm, 'ep2d_simpace_')[0]
        
        #remove the run number suffix, add to the list
        new_string = new_string.split('_')[0]
        mot_order.append(new_string)
    
    return mot_order
    

def get_date(dcm_dirs):
    """This function takes in a list of ordered dicom directories and returns the date
    of acquisition for the first in the list"""

    #load the first dicom from the first dir in the list
    test_dcm = sorted(glob(pjoin(dcm_dirs[0],'IM-*.dcm')))[0]
    hdr = dicom.read_file(test_dcm)
    date = hdr.AcquisitionDate

    return date


def rename_niftis(prefix,nii_dir,numvols):
    """This function renames niftis in nii_dir with a new prefix instead of the
    nipype given prefix (from dicom header)"""

    #get the nifti files in the directory
    niftis = sorted(glob(pjoin(nii_dir,'*.nii')))

    vols = range(numvols+1,(len(niftis)+numvols+1)) #get the number of volumes total to create new files

    #loop through them and rename the prefix
    for fidx, fname in enumerate(niftis):
        
        new_fname = """%s-%04d.nii""" %(prefix,vols[fidx])
        # print new_fname
        
        #change the filenames
        shutil.move(fname,pjoin(nii_dir,new_fname))


def move_first_vols(new_dcmdir,init_dcmdir,numvols=4):
    """This function moves a given number of volumes from new_dcmdir to init_dcmdir"""

    for vol in range(1,numvols+1):
        dcms = glob(pjoin(new_dcmdir,'IM-*-%04d.dcm' %(vol)))
        dcm = dcms[0]
        #TEST that dcm only has one value
        if not len(dcms) == 1:
            raise ValueError("more than on dicom match for moving")

        mv_command = """mv %s %s""" %(dcm,init_dcmdir)
        os.system(mv_command)

#-----------------------------------------------------------------------------
# Main Script
#-----------------------------------------------------------------------------
DESPO_SIMPACE_DIR = '/home/despo/simpace/subject_1_data/'
SUB_NUM = 1
NB_DISCARD_VOL = 4 #number of initial volumes to discard

def main(argv = sys.argv):
    """
    argv[1]: the top subject directory. 
             optional, if not given defaults to DESPO_SIMPACE_DIR
    """
    #get list of dicoms
    if len(argv) >= 2:
        # first argument is the top subject directory
        sub_dir = argv[1]
    else:
        sub_dir = DESPO_SIMPACE_DIR 
    assert check_dir(sub_dir, ['exists'])

    numvols = NB_DISCARD_VOL 
    sub = SUB_NUM
    
    sess_dirs = glob(pjoin(sub_dir,'ImageData*'))
    sess_dirs = sort_ctime(sess_dirs)

    #set up dicom conversion for later
    matlab.MatlabCommand.set_default_paths('/usr/local/matlab-tools/spm/spm8')
    dicom_convert = spmu.DicomImport()

    for sessidx, sess in enumerate(sess_dirs):

        print "working on session " + sess
        #setup the output directory for the new file names,
        out_dir = pjoin(sub_dir,'rename_files','sess%02d' %(sessidx+1))
        
        #rm the directory and contents if it already exists
        # !!!!!!! DEBUG MODE : remove dir if exists 
        if check_dir(out_dir, ['exists']): # equivalent to osp.isdir(out_dir)
            try:
                shutil.rmtree(out_dir)
            except:
                print "cannot remove " + out_dir

        #only run if this directory does not exist yet
        if not check_dir(out_dir, ['exists']):
            try:
                os.makedirs(out_dir, 0o770) 
            except:
                print "cannot make " + out_dir

        #get the two files for each motion type
        none_dirs = glob(pjoin(sess,'*/Danzone_Testing*/ep2d_simpace_NONE*'))
        low_dirs = glob(pjoin(sess,'*/Danzone_Testing*/ep2d_simpace_LOW*'))
        med_dirs = glob(pjoin(sess,'*/Danzone_Testing*/ep2d_simpace_MED*'))
        high_dirs = glob(pjoin(sess,'*/Danzone_Testing*/ep2d_simpace_HIGH*'))

        #select the file without prescan normalization
        fnone = get_nonorm(none_dirs)[0]
        flow = get_nonorm(low_dirs)[0]
        fmed = get_nonorm(med_dirs)[0]
        fhigh = get_nonorm(high_dirs)[0]

        #put all 4 runs into a list for sorting
        runs = [fnone, flow, fmed, fhigh]
        runs_sort = sort_acqtime(runs)

        #copy dicoms to have new directory
        for runidx, runname in enumerate(np.array(runs_sort)):

            #create the new dir name, and make it
            new_fname = 'sub%02d_sess%02d_run%02d' %(int(sub),sessidx+1,runidx+1)
            new_runname = pjoin(out_dir,new_fname)
            new_dcmdir = pjoin(new_runname,'dicoms')
            print " \t working on run : " + runname # debug
            print " \t run name : " + new_runname # debug
            print " \t new_dcmdir : " + new_dcmdir

            #-------- Directory not created yet - cannot check 
            #- if not (isuser_writeable(new_runname) and isuser_executeable(new_runname)):
            #-     # attempt to change permissions
            #-     os.chmod(new_runname, 0o770)
            #- if not (isuser_writeable(new_dcmdir) and isuser_executeable(new_dcmdir)):
            #-     # attempt to change permissions
            #-     os.chmod(new_dcmdir, 0o770)

            #copy dir with the new name
            try:
                # which permissions are given in the copy ?
                shutil.copytree(runname,new_dcmdir)
            except:
                print "cannot copytree " + runname + " in " + new_dcmdir

            init_dcmdir = pjoin(new_dcmdir,'first_vols')
            if not os.path.exists(init_dcmdir):
                try:
                    os.makedirs(init_dcmdir, 0o770)
                except:
                    st = os.stat(init_dcmdir)
                    print "cannot create 'first_vols' in " + new_dcmdir + str(st.st_mode)
                    raise

            #move the first four dcms to a new dir
            move_first_vols(new_dcmdir,init_dcmdir,numvols=numvols)
            
            #do dicom conversion
            #first create a nifti directory
            nii_dir = pjoin(new_runname,'data')
            if not check_dir(nii_dir):
                os.makedirs(nii_dir, 0o770)
            else:
                print "nifti directory " + nii_dir + " already exists"
                raise

            dcm_files = sorted(glob(pjoin(new_dcmdir,'*.dcm')))
            dicom_convert.inputs.in_files = dcm_files
            dicom_convert.inputs.output_dir = nii_dir
            dicom_convert.run()
            
            #rename the niftis in this directory 
            rename_niftis(new_fname,nii_dir,numvols)
            

        #get info for the params file
        #date of data collection (nb: this assumes all runs were collected in the same day)
        scan_date = get_date(runs_sort)

        #order of motion conditions
        mot_order = get_motorder(runs_sort)
        
        #save params file
        json_fname = pjoin(out_dir,'sub%02d_sess%02d_params' %(int(sub),sessidx+1))
        params = dict(date = scan_date,
                      motion = mot_order)
        with open(json_fname, 'w') as outfile:
            json.dump(params, outfile)
        #1/0


if __name__ == '__main__':
    main()
