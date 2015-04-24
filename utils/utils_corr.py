from __future__ import print_function
import nibabel as nib
import numpy as np
import os
import os.path as osp
import matplotlib.pyplot as plt
import nipy
import json
import scipy.linalg as lin
from warnings import warn

from numpy.testing import (assert_array_almost_equal, assert_almost_equal, assert_array_equal, assert_equal)



# assume both in MNI space - resample 
from nipy.algorithms.resample import resample_img2img, resample
import nipy.core.api as capi # vox2mni
from nipy.core.image import image
from nipy.core.reference import coordinate_map as cmap
from nipy.core.reference.coordinate_system import CoordSysMaker


def extract_signals(rois_dir, roi_prefix, fn_img4d, mask=None):
    """
    Extract signals from list of 4d nifti
    
    inputs
    ------
    rois_dir: str
    roi_prefix: str
    fn_img4d: str
    mask: str or nibabel image
    """
    
    IMGDIM = 3
    TINY = np.finfo('float').eps * 1000
    
    img4d = nib.load(fn_img4d)
    aff_img = ucr.xyz_affine(img4d.get_affine(), xyz=[0,1,2])
    
    # 1- Check roi affine are all the same
    aal_files = gb.glob(osp.join(rois_dir, roi_prefix))
    roi_imgs = {}
    for fn in aal_files:
        roi_imgs[fn] = nib.load(fn)
    roi_affs = np.asarray([img.get_affine()[:IMGDIM,:IMGDIM] for img in roi_imgs.values()])
    aff_roi = roi_affs[0].copy()
    roi_affs -= aff_roi
    norms_affs = [lin.norm(roi_affs[i], 'fro') for i in range(len(roi_affs))]
    assert lin.norm(np.asarray(norms_affs))  < TINY 
        
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
        elif hasattr(mask, get_data):
            data_mask = img_mask.get_data()
        else:
            data_mask = np.asarray(mask)
        data_mask = data_mask.astype('bool')
    
    
    # 3- extract and return signals
    signals = {}
    issues = {}
    for roi_name, roi in roi_imgs.iteritems():
        print(roi_name)
        ijk_roi = np.asarray(np.where(roi.get_data().astype('bool')))
    
        if translate:
            ijk_roi += translation
        
        print(ijk_roi.min(axis=1), ijk_roi.max(axis=1))
        ## is this translation ok with the shape of the image to sample from?
        img3dshape = np.asarray(img4d.shape)[:IMGDIM].reshape(IMGDIM,1)
        print(img3dshape)
        
        if not (np.all(ijk_roi >= 0) and np.all((img3dshape - ijk_roi) > 0)):
            print("could not sample roi :", roi_name)
            signals[roi_name] = None
            issues[roi_name] = 'ijk_roi out of image shape'
            continue
            
        # see if all ijk are within the mask?
        if mask != None:
            # not all voxels in the mask ? 
            vox_in_mask = np.all(data_mask[ijk_roi], axis=0)
            if not np.all(data_mask[ijk_roi]):
                issues[roi_name] = ('not all voxels in mask', len(vox_in_mask==True))
                print("some voxels not in mask")
        else:
            vox_in_mask = np.ones(ijk_roi.shape[1]).astype('bool')
        
        signals[roi_name] = (img4d.get_data()[ijk_roi[0][vox_in_mask], 
                                              ijk_roi[1][vox_in_mask], 
                                              ijk_roi[2][vox_in_mask], :]).mean(axis=0)
        
    return signals


def xyz_affine(big_aff, xyz=[0,1,2], debug=0):
    """
    take a big affine, 
    and returns the affine constructed with the indices xyz
    future : could extend to make it return a bigger affine than input aff (big_aff)
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



from nipy.labs.mask import compute_mask #, compute_mask_files

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


