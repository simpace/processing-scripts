from __future__ import print_function
import nibabel as nib
import numpy as np
import os.path as osp
import glob as gb
import scipy.linalg as lin
from six import string_types
from collections import OrderedDict

# import matplotlib.pyplot as plt
# import nipy
# import json
# from warnings import warn
# from numpy.testing import (assert_array_almost_equal, assert_almost_equal, \
#                               assert_array_equal, assert_equal)

#---  assume both in MNI space - resample 
# import nipy.core.api as capi # vox2mni
# from nipy.core.image import image
from nipy.algorithms.resample import resample # ,resample_img2img 
from nipy.core.reference import coordinate_map as cmap
from nipy.core.reference.coordinate_system import CoordSysMaker
from nipy.labs.mask import compute_mask #, compute_mask_files

TINY = np.finfo('float').eps * 1000

def project_filter(data, bandf, dt):
    """
    bandf : list or tuple
        [lowfreq, highfreq]
    """
    nt = data.shape[0]
    lowf = bandf[0]
    higf = bandf[1]
    tim = np.linspace(0, (nt-1)*dt, nt)
    #cosBF = dm._cosine_drift(1./lowf, frametimes)
    #cosBF = np.loadtxt('./dctmtx_114_114.txt')
    cosBF = _cosine_drift((2*dt), tim)
    #print(cosBF.shape)
    
    #order_h = int(np.floor(2*T) * higf)
    #order_l = int(np.floor(2*T) * lowf)
    order_h = int(np.floor(2*nt*higf*dt))
    order_l = int(np.floor(2*nt*lowf*dt))
    print(order_h, order_l)
    
    if order_h == 0:
        XBF = cosBF[:,:order_l]
    elif order_l == 0:
        XBF = cosBF[:,order_h:]
    else:
        XBF = np.hstack((cosBF[:,:order_l], cosBF[:,order_h:]))
    
    return R_proj(XBF, data)


def _create_bandpass_bf(npoint, dt, bandf, exclude=False):
    """
    npoint: int
        number of points in data to be filtered (ie. time dimension)
    bandf : list or tuple
        [lowfreq, highfreq]
    dt: number 
        the sampling time in seconds, eg. 2.4
    exclude: bool
        True will return basis functions between lowf and highf
        False will return bf between [0,lowf] and [higf,1/dt]
    """
    
    nt = npoint
    lowf = bandf[0]
    higf = bandf[1]
    tim = np.linspace(0, (nt-1)*dt, nt)
    cosBF = _cosine_drift((2*dt), tim)
    
    order_h = int(np.floor(2*nt*higf*dt))
    order_l = int(np.floor(2*nt*lowf*dt))
    # for debug: 
    # print(order_h, order_l)

    if lowf >= higf:
        raise(ValueError, (lowf, higf))

    # check at least one order positive
    if (order_h == 0) and (order_l == 0):
        raise(ValueError, (order_h, order_l))


    if order_h == 0:
        XBF = cosBF[:,:order_l]
    elif order_l == 0:
        XBF = cosBF[:,order_h:]
    else:
        XBF = np.hstack((cosBF[:,:order_l], cosBF[:,order_h:]))
    
    return XBF, (order_h, order_l)



#- in the future : replace this by import the right version of nipy -#
def _cosine_drift(period_cut, frametimes):
    """Create a cosine drift matrix with periods greater or equals to period_cut

    Parameters
    ----------
    period_cut: float 
         Cut period of the low-pass filter (in sec)
    frametimes: array of shape(nscans)
         The sampling times (in sec)

    Returns
    -------
    cdrift:  array of shape(n_scans, n_drifts)
             cosin drifts plus a constant regressor at cdrift[:,0]

    Ref: http://en.wikipedia.org/wiki/Discrete_cosine_transform DCT-II
    """
    len_tim = len(frametimes)
    n_times = np.arange(len_tim)
    hfcut = 1./ period_cut # input parameter is the period  

    dt = frametimes[1] - frametimes[0] # frametimes.max() should be (len_tim-1)*dt    
    order = int(np.floor(2*len_tim*hfcut*dt)) # s.t. hfcut = 1/(2*dt) yields len_tim
    cdrift = np.zeros((len_tim, order))
    nfct = np.sqrt(2.0/len_tim)
    
    for k in range(1, order):
        cdrift[:,k-1] = nfct * np.cos((np.pi/len_tim)*(n_times + .5)*k)
    
    cdrift[:,order-1] = 1. # or 1./sqrt(len_tim) to normalize
    return cdrift


#- use nilearn to get to get the rois and check all affines are the same

#- get  4d image - return array? iterator over nibabel images? 

#-----------------------------------------------------------------------------#
#-  This function should be a wrapper over nilearn stuff -#
#-----------------------------------------------------------------------------#


def _check_images_compatible(roi_files, nifti_image, imgdim=3):

    if isinstance(roi_files, string_types):
        roi_files = [roi_files]

    nb_rois = len(roi_files)

    # construct a list of affine, check the dimension of the list 
    roi_affs =  np.asarray([img.get_affine() for img in 
                                        [nib.load(fn) for fn in roi_files] ])

    assert roi_affs.shape == (nb_rois, imgdim+1, imgdim+1), \
            "roi_affs.shape unexpected: {}".format(roi_affs.shape)

    # check all the same affine as the first one
    aff_roi = roi_affs[0].copy()
    roi_affs -= aff_roi
    norms_affs = [lin.norm(roi_affs[i], 'fro') for i in range(len(roi_affs))]
    assert lin.norm(np.asarray(norms_affs))  < TINY, "roi affines not all the sames"

    # now check that the image has the same affine
    # check we can get img affine - THIS ASSUMES A 3D AFFINE - check to be done


    if isinstance(nifti_image, string_types):
        try: 
            nifti_image = nib.load(nifti_image)
        except:
            raise ValueError, "niti image, are you sure ? {}".format(nifti_image)

    try:
        img_aff = nifti_image.get_affine()
    except:
        raise ValueError, "check nifti_image is an nifti image"

    assert np.allclose(aff_roi, img_aff), \
           "aff_roi : {} and img_aff: {} not close enough".format(aff_roi, img_aff)

    return True


def _get_and_sort_roi_files(rois_dir, roi_prefix):

    roi_files = gb.glob(osp.join(rois_dir, roi_prefix))
    #print(roi_files)

    # fails if roi_files == []
    assert roi_files, " roi_files empty {}".format(roi_files)
    # sort the files to have some consistency in the behaviour
    roi_files.sort()

    return roi_files


def _yield_roi_mask(roi_files, imgdim=3, roi_threshold=0.5):
    """
    this works even if there are some NaNs in the roi because
    NaNs compare to float yield False : the roi will not
    contain NaNs
    """

    for fn_roi in roi_files:
        label = osp.splitext(osp.basename(fn_roi))[0]
        img = nib.load(fn_roi)
        roi_mask = np.asarray(img.get_data().astype(float) > roi_threshold)
        assert roi_mask.dtype == np.dtype('bool')

        yield label, roi_mask

def extract_signals(rois_dir, roi_prefix, fn_img4d, mask=None, 
                                            minvox=1, verbose=False):
    """
    Extract signals from list of 4d nifti 
    
    Parameters
    -----------
    rois_dir: str
    roi_prefix: str
    fn_img4d: str or img
    mask: str or nibabel image
    minvox: minimum number of voxel in sampled region

    returns:
    --------
    signals: dict
        keys are the names of the regions (from filenames)
        values: None, or the signal in the region
    issues: dict
        keys are the names of the regions (from filenames)
        values: None if there is no issue, a tuple otherwise
    """

    IMGDIM = 3
    roi_files = _get_and_sort_roi_files(rois_dir, roi_prefix)
    
    if verbose: print(roi_files[:3],'\n')


    if isinstance(mask, string_types):
        mask = nib.load(mask)
    if isinstance(fn_img4d, string_types):
        img4d = nib.load(fn_img4d)
    else:
        img4d = fn_img4d

    assert  _check_images_compatible(roi_files, img4d, imgdim=IMGDIM)


    if mask is not None:
        assert np.allclose(mask.get_affine(), img4d.get_affine()), \
                "mask {} and img4d {} affine not close ".format(
                      mask.get_affine(), img4d.get_affine())
        mask_arr = np.asarray(mask.get_data()).astype(bool)
        # some checks on mask ?
        nb_vox_mask = mask_arr.sum()
        nb_vox_img  = np.asarray(mask_arr.shape).prod()
        if verbose:
            print("In mask:{}, out:{}".format(nb_vox_mask, nb_vox_img))
        assert nb_vox_mask > .20*nb_vox_img

    img_arr = np.asarray(img4d.get_data())
    
    # initialize dicts 
    signals = {}
    info = {}
    issues = {}
    for label, roi_msk in _yield_roi_mask(roi_files):
        if verbose: print("nb voxels in mask of {} is {} \n".format(label, roi_msk.sum()))
        if mask is not None:
            this_msk = np.logical_and(roi_msk, mask_arr)
        else: 
            this_msk = roi_msk

        this_msk_nb_vox = this_msk.sum()
        
        if this_msk_nb_vox < minvox:
            issues[label] = "this_msk.sum() < minvox: {} < {}".format(
                                    this_msk_nb_vox, minvox)
        else:
            # dimension should be (# of voxels in ROI and mask, time)
            signals[label] = img_arr[this_msk].mean(axis=0)
            info[label] = this_msk_nb_vox

    return signals, issues, info


def _dict_signals_to_arr(dsig):
    """
    Take a signal dictionary and turn the not None values into numpy array
    Suppose signals is a dict generated by extract_signals 
    """
    _keys = dsig.keys()
    #_keys = tuple([ k for k in _keys if dsig[k] != None ])
    t = dsig[_keys[0]].shape[0]
    n = len(_keys)
    arr = np.zeros((t, n), dtype='float')
    for idx, k in enumerate(_keys):
        arr[:,idx] = dsig[k][0]

    return (arr, _keys)


def extract_mvt_run(mvtfile, run_idx, run_len):
    """
    """
    mvt = np.loadtxt(mvtfile)
    assert mvt.ndim == 2
    assert mvt.shape[1] == 6
    assert (mvt.shape[0] % run_len) == 0
    
    return mvt[run_idx*run_len:(run_idx+1)*run_len,:]


#=====================================================
# general script
#=====================================================

from nilearn import masking as msk
from nilearn._utils import concat_niimgs
import shutil
import os 

# from nilearn.image.image import _compute_mean

def do_one_sess(DDIR, idx_sess):
    """
    take a session directory and extract the roi signals
    """
    runs_pat = ['*sess{:02d}_run{:02d}-0*'.format(idx_sess, idx_run) 
                                                for idx_run in [1,2,3,4]]
    runs = [gb.glob(osp.join(DSRA, pat)) for pat in runs_pat]
    for run in runs: run.sort()

    # compute session wide mask
    #-----------------------------------------------------

    sess_mask = msk.compute_multi_epi_mask(runs, lower_cutoff=0.2, upper_cutoff=0.85,
                           connected=True, opening=2, threshold=0.5)
    #compute_epi_mask(runs[0], opening=1, connected=True)
    
    # check mask is reasonable - how ???

    for idx_run, run in enumerate(runs):
        info = do_one_run(run, idx_run, sess_mask)


def do_one_run(file_names, idx_run, mask):
    """
    """

    # extract base directory from first file
    ddir = osp.dirname(file_names[0])
    droi = osp.join(ddir, 'registered_files')
    dsig = osp.join(ddir, 'extracted_signals')
    fn_sig = osp.join(dsig, 'signal_run_{:02d}'.format(idx_run))

    # extract signals and save them in preproc/roi_signals
    #-----------------------------------------------------

    roi_prefix = 'rraal_*.nii'
    run_4d = concat_niimgs(file_names, ensure_ndim=4)
    signals, _issues, _info = extract_signals(droi, roi_prefix, run_4d, 
                                                        mask=mask, minvox=1)   
    arr_sig, labels = _dict_signals_to_arr(signals)

    # save in 
    if osp.isdir(dsig):
        try:
            shutil.rmtree(dsig)
        except:
            # try change permissions ?
            raise ValueError, 'cannot rm {}, check perm.'.format(dsig)
    else:
        try:
            os.makedirs(dsig, 0o770)
        except:
            raise ValueError, 'cannot create {}, check perm.'.format(dsig)

    np.savez(fn_sig, arr_sig, labels) 


    # construct matrix of counfounds
    #-----------------------------------------------------


    # filter and compute correlation
    #-----------------------------------------------------

    # save correlation




#----------------------------------------------------------------------------#
# A list of functions that return a bool array with True where signal ok
# all have inputs arr
#----------------------------------------------------------------------------#
def _var_not_zero(arr):
    """
    check that variance of signal is strong enough
    """
    return np.var(arr, axis=0) > TINY 


def _check_roi_signals(arr, kkeys, check_funcs = [_var_not_zero]):
    """
    logical and with all functions in check_funcs
    """
    kkeys = np.asarray(kkeys)
    goods = np.ones((len(kkeys),),dtype='bool',)
    for func in check_funcs:
        goods = np.logical_and(goods, func(arr))

    return arr[:,goods], kkeys[goods] 



def xyz_affine(big_aff, xyz=[0,1,2], debug=0):
    """
    take a big affine, 
    and returns the affine constructed with the indices xyz
    could be extended to make it return a bigger affine than input big_aff
    """
    
    tmp = big_aff[xyz,:][:,xyz]
    trans = big_aff[xyz,-1].reshape((tmp.shape[0],1))
    last_row = np.zeros((1,tmp.shape[1]+1))
    last_row[0,-1] = 1.
    tmp = np.hstack((tmp,trans))
    tmp =  np.vstack((tmp, last_row))
    if debug: print('\nold affine:\n',big_aff, '\nnew affine:\n', tmp)

    return tmp


def resamp(nipy_img, onto_aff, onto_shape, order=3, cval=-1, w2wmap=None, debug=0):
    """
    utility function to resample an image
    input: 
	nipy image
    	onto affine
	onto shape
	order : scipy resample interpolation order
	cval: value outside when resampling to bigger box
    """
    arraycoo = 'ijklmnopq'[:len(onto_shape)]
    spacecoo = 'xyztrsuvw'[:len(onto_shape)] 
    if debug: print('\narraycoo: ', arraycoo,  
                     '\nspacecoo: ', spacecoo, '\nonto_aff\n', onto_aff)
    
    dmaker = CoordSysMaker(arraycoo, 'generic-array')
    rmaker = CoordSysMaker(spacecoo, 'generic-scanner')
    cm_maker = cmap.CoordMapMaker(dmaker, rmaker)
    cmap_out = cm_maker.make_affine(onto_aff)
    if debug: print('cmap_out:\n',cmap_out)
    
    if w2wmap == None: w2wmap = np.eye(onto_aff.shape[0])
    if debug: print('w2wmap:\n',w2wmap)
        
    return resample(nipy_img, cmap_out, w2wmap, onto_shape, order=order, cval=cval)
    
def ioflab(l,labels):
    return [ii for (ii,lab) in enumerate(labels) if (l in lab[0])]


def inter_run_mask(fnames, hist_m=.3, hist_M=.8):
    """
    create the inter runs mask : 
        - first create a mask of all constant voxels on axis 3
        - second create a mask of the mean 
        - 
    """
    EPSFCT = 10000. #  factor: multiply numpy epsilon by it
    imgs = [nib.load(fn) for fn in fnames]
    # initialize mask to ones
    mask = np.ones(imgs[0].shape[:3])
    for img in imgs:
        arrstd = np.asarray(img.get_data().std(axis=3))
        mask = np.logical_and(mask, arrstd > EPSFCT*np.finfo(arrstd.dtype).eps)
        mask_mean  = compute_mask(img.get_data().mean(axis=3), None, \
                            m=hist_m, M=hist_M, cc=True, opening=2, exclude_zeros=False)
        mask = np.logical_and(mask, mask_mean)
        
    return mask
        
# ----- Testing--------------------------------------- 
#arrstd = np.asarray(fimg.get_data().std(axis=3))
#arrstd.mean(), arrstd.min(), arrstd.max()
#print arrstd.dtype
#np.finfo(arrstd.dtype).eps
#mask = arrstd > 10*np.finfo(arrstd.dtype).eps
#print mask.dtype
#mask.sum(), np.prod(mask.shape)
# ----- Testing-------------------------------------- 

from nipy.modalities.fmri import design_matrix as dsgn_mtrx

SKIPTR=0
SK4RUNLENGH=196

def make_design_mtx(spm_mvt_file, hfcut=128, skip_TR=SKIPTR, normalize=True, to_add=None):
    """
    create the design matrix
    """
    mvt_arr = np.loadtxt(spm_mvt_file)
    mvt_arr = mvt_arr[skip_TR:,:]
    mvt_lab = ['tx','ty','tz','rx', 'ry', 'rz'] 
    
    assert mvt_arr.shape == (SK4RUNLENGH,6)
    if to_add is not None:
        assert to_add.shape[0] == SK4RUNLENGH
        addlab = [ 'a'+str(idx) for idx in range(to_add.shape[1])]
        mvt_arr = np.hstack((mvt_arr, to_add))
        mvt_lab = mvt_lab + addlab
    

    nbframes = mvt_arr.shape[0]
    frametimes = np.linspace(0.0, (nbframes-1)*2.0, num=nbframes)
    
    X, labels = dsgn_mtrx.dmtx_light(frametimes, paradigm=None, hrf_model='canonical',
              drift_model='cosine', hfcut=hfcut, 
              drift_order=1,
              fir_delays=[0], add_regs=mvt_arr, add_reg_names=mvt_lab)

    if normalize:
    	X -= X.mean(axis=0)
	Xstd = X.std(axis=0)
	non_zero_col = Xstd > np.finfo(float).eps
	zero_col = np.logical_not(non_zero_col)
	X = np.hstack((X[:,non_zero_col]/Xstd[non_zero_col], X[:, zero_col])) 
    
    return X, labels


#----------------------------------------------
# Projection functions

def proj(X, Y=None):
    # project onto the space of X
    #print matrix_rank(X)
    [u, s, vt] = lin.svd(X, full_matrices=False)
    tol = s.max() * max(X.shape) * np.finfo(s.dtype).eps
    nz = np.where(s > tol)[0]
    #print nz
    if Y is None:
        return u[:,nz].dot(u[:,nz].T)
    else:
        return u[:,nz].dot(u[:,nz].T.dot(Y))
    
def R_proj(X, Y=None):
    # project onto the space orthogonal to X (RY) 
    if Y is None:
        return np.eye(X.shape[0]) - proj(X, Y) 
    else:
        return Y - proj(X, Y)

#----------------------------------------------

def process_one_roi(roi_idx, fmri_data, roi_mask, X, skip_TR=SKIPTR):
    """
    return the signal in the ROI indexed by roi_idx
    after regressing out X
    """
    # mask_roi 
    # roi_mask = np.logical_and((atlas_data == roi_idx), (fmri_mask==1))
    
    Nvox = roi_mask.sum()    
    roi_data = fmri_data[roi_mask, :].T # .T : make it time x voxels
    roi_data = roi_data[skip_TR:,:]
    #print "roi_data.shape " , roi_data.shape
    roi_data = R_proj(X, roi_data)
    return roi_data.mean(axis=1), Nvox
    
    
def mask_of_roi(roi_idx, fmri_mask, atlas_data):
    """
    return the mask voxels in roi     """
    idx_vox = np.logical_and((atlas_data == roi_idx), fmri_mask)
    return idx_vox

def extract_roi_signals(atlas, fmri, labels, roiMinSize=0):
    """

    """
    NbLabels = len(labels)
    extracted = np.zeros((fmri.shape[-1], NbLabels))
    extracted_labels = []
    large_sizes = []
    too_small = []
    small_sizes = []

    out_idx = 0

    for ii,lab in enumerate(labels):
        assert lab[0].shape[0] == 1
        idx = int(lab[0]) # one value 
        roi = (atlas == idx)
        NbVox = roi.sum()

        if NbVox > roiMinSize:
            extracted_labels.append(lab)
            large_sizes.append(NbVox)
            signal = fmri[roi,:].mean(axis=0)
            extracted[:, out_idx] = signal
            out_idx += 1
        else:
            too_small.append(lab)
            small_sizes.append(NbVox)
           
    goods = zip(extracted_labels, large_sizes)
    bads = zip(too_small, small_sizes)
    return extracted[:,:out_idx], goods, bads


def lab_names(labels):
    return [lab[1][1] for lab in labels]

def lab_idx(labels):
    return [int(lab[0]) for lab in labels]

def glab_enum(labels):
    current = 0
    while current < len(labels):
        yield lab_idx([labels[current]])[0], lab_names([labels[current]])[0]
        current += 1
	
def lab_enum(labels):
    return [(lab_idx([labels[current]])[0], lab_names([labels[current]])[0]) 
            for current in range(len(labels))]
	

# correlations = []
# 
# idx_1 = 2
# idx_2 = 3
# 
# for run_idx in range(EPIRUNSNB):
#     X, _ = make_design_mtx(spm_mvt_files[run_idx], hfcut=128, skip_TR=SKIPTR)
#     roi_mask = mask_of_roi(idx_1, mask_all_runs, dratlas)
#     sig_1, _ = process_one_roi(idx_1, np.asarray(fmri_runs[run_idx].get_data()), \
#                                   roi_mask, X, skip_TR=SKIPTR)
#     roi_mask = mask_of_roi(idx_2, mask_all_runs, dratlas)
#     sig_2, _ = process_one_roi(idx_2, np.asarray(fmri_runs[run_idx].get_data()), \
#                                   roi_mask, X, skip_TR=SKIPTR)
#     
#     if (sig_1 is not None) and (sig_2 is not None):
#         correlations.append(np.corrcoef(sig_1, sig_2)[0,1])

#def make_all_roi_roi_corr(fruns, mvts, mask_all_runs, dratlas, labels, hfcut=128, skip_TR=SKIPTR):
#    """
#    """
#    MINROISIZE = 5
#    nb_labels = len(labels)
#    nb_runs = len(fruns)
#    results_corr = np.ones((nb_labels, nb_labels, nb_runs)) + 10. # put ten : not a correlation
#    
#    Xs = [make_design_mtx(mvt, hfcut=hfcut, skip_TR=SKIPTR)[0] for mvt in mvts]
#    dfruns = [np.asarray(img.get_data()) for img in fruns]
#    resdic = {}
#    bad_rois = []
#    
#    for idx1,lab1 in enumerate(labels):
#        
#        region_name1 = lab1[1][1]
#        if region_name1 in bad_rois:
#            continue
#        
#        mask_roi_1 = mask_of_roi(idx1, mask_all_runs, dratlas)
#        nb_vox_ROI_1 = mask_roi_1.sum()
#        
#        if nb_vox_ROI_1 < MINROISIZE:
#            bad_rois.append(region_name1)
#            continue
#        else:
#            print idx1, region_name1, nb_vox_ROI_1, ": ... " 
#            
#            for idx2,lab2 in enumerate(labels[idx1:]):
#                # compute correlations for idx1, idx2
#                region_name2 = lab2[1][1]
#                if region_name2 in bad_rois:
#                    continue
#
#                mask_roi_2 = mask_of_roi(idx2, mask_all_runs, dratlas)
#                nb_vox_ROI_2 = mask_roi_2.sum()
#                
#                if nb_vox_ROI_2 < MINROISIZE:
#                    bad_rois.append(region_name2)
#                    continue
#                else: 
#                    # we found a good pair 
#                    resdic[region_name1] = {}
#                    resdic[region_name1][region_name2] = [idx1,idx2]
#
#                    for run in range(nb_runs):
#                        sig_1, _ = process_one_roi(idx1, dfruns[run], mask_roi_1, \
#                                                   Xs[run], skip_TR=skip_TR)
#                        sig_2, _ = process_one_roi(idx2, dfruns[run], mask_roi_2, \
#                                                   Xs[run], skip_TR=skip_TR)
#                        results_corr[idx1, idx2, run] = np.corrcoef(sig_1, sig_2)[0,1]
#
#    
#    return results_corr, resdic, bad_rois
#                
#        

def _extract_signals(rois_dir, roi_prefix, fn_img4d, mask=None, minvox=1):
    """
    Extract signals from list of 4d nifti 
    
    Parameters
    -----------
    rois_dir: str
    roi_prefix: str
    fn_img4d: str
    mask: str or nibabel image
    minvox: minimum number of voxel in sampled region

    returns:
    --------
    signals: dict
        keys are the names of the regions (from filenames)
        values: None, or the signal in the region
    issues: dict
        keys are the names of the regions (from filenames)
        values: None if there is no issue, a tuple otherwise
    """
    
    IMGDIM = 3
    
    if isinstance(fn_img4d, string_types):
        img4d = nib.load(fn_img4d)
    elif hasattr(fn_img4d, 'get_affine') and hasattr(fn_img4d, 'get_data'):
        img4d = fn_img4d 
    else:
        raise ValueError('I dont think {} can be thought as a 4d image',fn_img4d)

    aff_img = xyz_affine(img4d.get_affine(), xyz=[0,1,2])
    
    # 1- Check roi affine are all the same
    roi_files = gb.glob(osp.join(rois_dir, roi_prefix))
    nb_rois = len(roi_files)

    # make sure list is not empty
    assert roi_files

    roi_imgs = {}
    for fn in roi_files:
        roi_imgs[fn] = nib.load(fn)

    roi_affs = np.asarray([img.get_affine()[:IMGDIM,:IMGDIM] 
            for img in roi_imgs.values()])
    assert roi_affs.shape == (nb_rois, IMGDIM, IMGDIM), \
            "roi_affs.shape unexpected: {}".format(roi_affs.shape)
    
    aff_roi = roi_affs[0].copy()
    roi_affs -= aff_roi
    norms_affs = [lin.norm(roi_affs[i], 'fro') for i in range(len(roi_affs))]
    assert lin.norm(np.asarray(norms_affs))  < TINY, "roi affines not all the sames"

    # 2- check that rois in roi_dir are compatible with img4d
    # print(aff_roi[:IMGDIM,:IMGDIM], aff_img[:IMGDIM,:IMGDIM])
    assert lin.norm(aff_roi[:IMGDIM,:IMGDIM] - aff_img[:IMGDIM,:IMGDIM], 'fro')  < TINY
    # print(aff_roi, aff_img))
    
    # translate ?
    translate = False
    print(aff_img)
    full_aff_roi = roi_imgs[roi_imgs.keys()[0]].get_affine()
    print(full_aff_roi)
    roi2img = lin.inv(aff_img).dot(full_aff_roi)
    translation = (roi2img[:IMGDIM, -1]).reshape(IMGDIM, 1)
    print("translation: ",translation)
    if lin.norm(translation, 'fro') > TINY:
        translate = True
    
    # need to mask? get the mask - can be filename, nibabel image or array
    if mask != None:
        if isinstance(mask, string_types):
            data_mask = nib.load(mask).get_data()
        elif hasattr(mask, 'get_data'):
            data_mask = mask.get_data()
        else:
            data_mask = np.asarray(mask)
        data_mask = data_mask.astype('bool')
    
    
    # 3- extract and return signals
    signals = {}
    issues = {}
    for roi_name, roi in roi_imgs.iteritems():
        issues[roi_name] = False
        ijk_roi = np.asarray(np.where(roi.get_data().astype('bool')))
        print(roi_name)
        print(ijk_roi.shape)
    
        if translate:
            ijk_roi += translation
        
        ## is this translation ok with the shape of the image to sample from?
        img3dshape = np.asarray(img4d.shape)[:IMGDIM].reshape(IMGDIM,1)
        # print(img3dshape)
        
        if not (np.all(ijk_roi >= 0) and np.all((img3dshape - ijk_roi) > 0)):
            print("could not sample all roi :", roi_name)
            print(ijk_roi.min(axis=1), ijk_roi.max(axis=1))
            signals[roi_name] = None
            issues[roi_name] = 'ijk_roi out of image shape'
            continue
            
        # see if all ijk are within the mask? if not, sample region on mask
        if mask != None:
            # not all voxels in the mask ? 
            vox_in_mask = np.all(data_mask[ijk_roi], axis=0)
            print('vox_in_mask dim {}'.format(vox_in_mask.shape))
            if not np.all(data_mask[ijk_roi]):
                issues[roi_name] = ('not all voxels in mask', vox_in_mask.sum())
                print("some voxels not in mask")
        else:
            vox_in_mask = np.ones(ijk_roi.shape[1]).astype('bool')

        nvox = vox_in_mask.sum() # number of True == number of voxels in mask
        if nvox < minvox:
            signals[roi_name] = None
            issues[roi_name] = ('less than minvox = {:d} in roi'.format(minvox), nvox)
            continue
        
        signals[roi_name] = ((img4d.get_data()[ijk_roi[0][vox_in_mask], 
                                              ijk_roi[1][vox_in_mask], 
                                              ijk_roi[2][vox_in_mask], :]).mean(axis=0), nvox)
        
    return signals, issues
        


