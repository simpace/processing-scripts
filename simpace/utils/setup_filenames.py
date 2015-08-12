"""Script to setup dicom directories to standard format
"""


import os, sys, stat
import numpy as np
import re
import json
from glob import glob
from os.path import join as pjoin
import os.path as osp
import dicom
import shutil
import tempfile
from six import string_types

# export MATLAB_INSTALLED=1 
MATLAB_INSTALLED = os.environ.get('MATLAB_INSTALLED')
#if MATLAB_INSTALLED is not None:
import nipype.interfaces.spm.utils as spmu
import nipype.interfaces.matlab as matlab
#


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
    checklist: string or list of string
        with all the check to be performed
        ['exists', 'iswritable']
    returns:
    -------
    bool 
        a boolean if all checklist conditions are true
    """
    if isinstance(checklist, string_types):
        checklist = [checklist]

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
        else:
            raise ValueError, ' {} not implemented '.format(chck)

    return checkbool

def sort_ctime(file_list):
    """This function takes in a list of files and returns a sorted list of them 
    by their creation time"""
    
    files_sorted = sorted(file_list, key=os.path.getctime)

    return files_sorted


def sort_acqdate(file_list):
    """This function takes in a list of files and returns a sorted list of them 
    by the aquisition time in the dicom header from the localizer"""

    #initialize a list for the dates files
    dates = []
    
    #loop through the directories
    for fileidx, fname in enumerate (file_list):

        #load the first dicom 
        test_dcm = sorted(glob(pjoin(fname,'*','Danzone_Testing*','localizer*','IM-*.dcm')))[0]
        hdr = dicom.read_file(test_dcm)
        date = hdr.AcquisitionDate
        
        dates.append(date)

    #sort the dates and apply this to the directories
    sort = np.array(dates).argsort()
    files_sorted = np.array(file_list)[sort]

    return files_sorted.tolist()


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
    #if more than one match, find the one with the full number of volumes
    if not len(no_psnorm) == 1:
        no_psnorm = test_nvols(no_psnorm)

    return no_psnorm


def test_nvols(dcm_dirs, nvols_full=200):
    """This function reads in a list of directories and returns the one with a full
    run of data (i.e., the full number of volumes expected)"""

    #initialize a list for directories with full data
    full_data = []

    #loop through the dicom directories
    for dcmidx, dcm in enumerate (dcm_dirs):

        nvols = len(glob(pjoin(dcm,'IM-*.dcm')))
        print dcm, nvols
        if nvols == nvols_full:
            full_data.append(dcm)

    #TEST that there is only one dir being returned
    if not len(full_data) == 1:
        raise ValueError("more than one run with full data for a condition!")

    return full_data


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
        
        #change the filenames
        shutil.move(fname,pjoin(nii_dir,new_fname))


def move_first_vols(new_dir,init_dir,numvols=4):
    """This function moves a given number of volumes from new_dir to init_dir"""

    for vol in range(1,numvols+1):
        niis = glob(pjoin(new_dir,'*-*-%05d-%06d-*.nii' %(vol,vol)))
        nii = niis[0]
        
        #TEST that dcm only has one value
        if not len(niis) == 1:
            raise ValueError("more than on nifti match for moving")

        mv_command = """mv %s %s""" %(nii,init_dir)
        os.system(mv_command)


def rm_and_create(rm_dir, perm=0o770):

    # 0o770 : rwxrws___  
    def remove_readonly(func, path, excinfo):
        os.chmod(path, stat.S_IWUSR)
        func(path)

    if os.path.isdir(rm_dir):
        try:
            # if fails, onerror tries to change permissions
            shutil.rmtree(rm_dir, onerror=remove_readonly)
        except:
            raise ValueError, 'cannot rm {}, check perm.'.format(rm_dir)
    try:
        os.makedirs(rm_dir, perm)
    except:
        raise ValueError, 'cannot create {}, check perm.'.format(rm_dir)


def _check_glob_res(res, ensure=None, files_only=True):
    """
    Will check the result of a glob: files must exist
    ensure == 1: will check there is only one file and return a string
    """

    if isinstance(res, string_types):
        res = [res]

    if isinstance(res, list):
        if ensure is not None:
            assert len(res) == ensure, "res: {} len(res) {} != ensure {}".format(
                                        res, len(res), ensure)
        if files_only:
            for f in res:
                assert osp.isfile(f), " {} is not a file".format(f)
    
    if ensure == 1:
        res = res[0]

    return res

        

#-----------------------------------------------------------------------------
# Main Script
#-----------------------------------------------------------------------------
DESPO_SIMPACE_DIR = '/home/despo/simpace/'
SUB_NUM = 1
NB_DISCARD_VOL = 4 #number of initial volumes to discard


def main(argv = sys.argv):
    """
    argv[1]: a string to specify a session, eg: sept_17, or all
    argv[2]: the top directory. 
             optional, if not given defaults to DESPO_SIMPACE_DIR
    argv[3]: the subject number.
             optional, if not given, defaults to SUB_NUM
    """
    if not MATLAB_INSTALLED:
        print("MATLAB_INSTALLED is {}:doing nothing".format(MATLAB_INSTALLED))
        #return None
    
    # first argument is the top subject directory
    if len(argv) >= 2:
        sess_specified = argv[1]
    else:
        sess_specified = ''
    if sess_specified == 'all':
        sess_specified = ''

    if len(argv) >= 3:
        # second argument is the top data directory
        data_dir = argv[2]
    else:
        data_dir = DESPO_SIMPACE_DIR 

    if len(argv) >= 4:
        # thrid argument is the subject directory
        sub_dir = pjoin(data_dir,'subject_%s_data' %(argv[3]))
    else:
        sub_dir = pjoin(data_dir,'subject_%s_data' %(SUB_NUM))
    
    assert check_dir(sub_dir, ['exists'])

    numvols = NB_DISCARD_VOL 
    sub = SUB_NUM

    #set up dicom conversion for later
    matlab.MatlabCommand.set_default_paths('/usr/local/matlab-tools/spm/spm8')
    dicom_convert = spmu.DicomImport()


    #CONVERT AND RENAME EPIs
    #get list of epi dicoms
    sess_dirs = glob(pjoin(sub_dir,'ImageData*'+sess_specified+'*'))
    sess_dirs = sort_acqdate(sess_dirs)

    for sessidx, sess in enumerate(sess_dirs,start=1):
        
        print "working on session " + sess
        #setup the output directory for the new file names,
        out_dir = pjoin(data_dir,'rename_files','sub%02d'%(sub),'sess%02d'%(sessidx))
        
        #rm the directory and contents if it already exists
        # !!!!!!! DEBUG MODE : remove dir if exists 
        if check_dir(out_dir, ['exists']): # equivalent to osp.isdir(out_dir)
            try:
                shutil.rmtree(out_dir)
            except:
                # attempt to change permissions ?
                # os.chmod(out_dir, 0o770)
                print "cannot remove : " + out_dir
                print "please run 'chmod -R 770 on : " + out_dir

        #only run if this directory does not exist yet
        if not check_dir(out_dir, ['exists']):
            try:
                os.makedirs(out_dir, 0o770) 
            except:
                print "cannot make " + out_dir
                raise
        else:
            print " directory : " + out_dir + " already exists"

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
        for runidx, runname in enumerate(np.array(runs_sort),start=1):

            #create the new dir name, and make it
            new_fname = 'sub%02d_sess%02d_run%02d' %(int(sub),sessidx,runidx)
            new_runname = pjoin(out_dir,new_fname)
            print " \t working on run : " + runname # debug
            print " \t run name : " + new_runname   # debug

            #get the list of dcm files for that run
            epidcm_files = sorted(glob(pjoin(runname,'*.dcm')))
            
            #create a nifti directory and do dicom conversion
            nii_dir = pjoin(new_runname,'niftis')
            if not check_dir(nii_dir):
                os.makedirs(nii_dir, 0o770)
            else:
                print "nifti directory " + nii_dir + " already exists"
                raise

            dicom_convert.inputs.in_files = epidcm_files
            dicom_convert.inputs.output_dir = nii_dir
            dicom_convert.run()

            #move the first niftis to a new dir
            init_nii_dir = pjoin(nii_dir,'first_vols')
            if not os.path.exists(init_nii_dir):
                try:
                    os.makedirs(init_nii_dir, 0o770)
                except:
                    st = os.stat(nii_dir)
                    print "cannot create 'first_vols' in " + nii_dir + \
                                                ' mode is :', str(st.st_mode)
                    raise
            move_first_vols(nii_dir, init_nii_dir, numvols=numvols)

            #rename the niftis 
            rename_niftis(new_fname,nii_dir,numvols)
            1/0

            #copy dir with the new name
            try:
                # which permissions are given in the copy ?
                shutil.copytree(runname, new_dcmdir)
            except:
                print "cannot copytree " + runname + " in " + new_dcmdir
                raise

            # The copy may have worked, but wrong permission passed. 
            # Try to 770 these.
            try:
                # walk from one level up the new_dcmdir
                for root, dirs, files in os.walk(new_runname):
                    for d in dirs:
                        os.chmod(pjoin(root,d),  0o770)
                    for f in files:
                        os.chmod(pjoin(root,f),  0o770)
            except:
                print "could not change permissions of " + new_dcmdir
                raise

            init_dcmdir = pjoin(new_dcmdir,'first_vols')
            if not os.path.exists(init_dcmdir):
                try:
                    os.makedirs(init_dcmdir, 0o770)
                except:
                    st = os.stat(new_dcmdir)
                    print "cannot create 'first_vols' in " + new_dcmdir + \
                                                ' mode is :', str(st.st_mode)
                    raise
            
            
            
        #get info for the params file
        #date of data collection (nb: this assumes all runs were collected in the same day)
        scan_date = get_date(runs_sort)

        #order of motion conditions
        mot_order = get_motorder(runs_sort)
        
        #save params file
        json_fname = pjoin(out_dir,'sub%02d_sess%02d_params.json' %(int(sub),sessidx))
        params = dict(date = scan_date,
                      motion = mot_order)
        with open(json_fname, 'wb') as outfile:
            json.dump(params, outfile)

        try:
            os.chmod(json_fname,  0o770)
        except:
            print "Could not chmod " + json_fname + " to 0o770 "
            raise


    #CONVERT AND RENAME ANAT
    #get list of epi dicoms
    anat_dir = glob(pjoin(sub_dir,'anatomical*','T1*'))
    #setup the output directory for anat
    out_adir = pjoin(data_dir,'rename_files','sub%02d'%(int(sub)),'anatomical')

    #rm the directory and contents if it already exists
    # !!!!!!! DEBUG MODE : remove dir if exists 
    if check_dir(out_adir, ['exists']): # equivalent to osp.isdir(out_adir)
        try:
            shutil.rmtree(out_adir)
        except:
            # attempt to change permissions ?
            # os.chmod(out_adir, 0o770)
            print "cannot remove : " + out_adir
            print "please run 'chmod -R 770 on : " + out_adir

    #only run if this directory does not exist yet
    if not check_dir(out_adir, ['exists']):
        try:
            os.makedirs(out_adir, 0o770) 
        except:
            print "cannot make " + out_adir
            raise
    else:
        print " directory : " + out_adir + " already exists"
    
    #create the new adir name for niftis
    nii_adir = pjoin(out_adir,'niftis')

    if not check_dir(nii_adir):
        os.makedirs(nii_adir, 0o770)
    else:
        print "anatomical nifti directory " + nii_adir + " already exists"
        raise

    anatdcm_files = sorted(glob(pjoin(anat_dir[0],'*.dcm')))
    dicom_convert.inputs.in_files = anatdcm_files
    dicom_convert.inputs.output_dir = nii_adir
    dicom_convert.run()

    #rename anat nifti
    old_anat_nii = glob(pjoin(nii_adir,'*.nii'))[0]
    new_anat_nii = pjoin(nii_adir,'sub%02d_anatomical.nii' %(int(sub)))
    shutil.move(old_anat_nii,new_anat_nii)
    

if __name__ == '__main__':
    main()
