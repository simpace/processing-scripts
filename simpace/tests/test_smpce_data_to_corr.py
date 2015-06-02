from __future__ import print_function, division
from nose.tools import raises

import scipy.io as sio
from numpy.testing import assert_allclose #, assert_array_equal
import os.path as osp
import glob as gb
from ..smpce_data_to_corr import get_params, process_all
from ..utils import _utils as ucr
from ..utils import setup_filenames as suf
from nilearn._utils import concat_niimgs

BASEDIR_JB = '/home/jb/data/simpace/data/rename_files'
BASEDIR_NX = '/home/despo/simpace/subject_1_data/rename_files'


@raises(IOError)
def test_process_all_assert():
    process_all('1')

def _get_matfile_data():

    # compare to Arielle's in the .mat file
    mat_file = "sub01_session_1_raw_ROI_timeseries.mat"
    mat_file = osp.join(osp.dirname(osp.realpath(__file__)), mat_file)
    # to_check.keys() == ['all_rois', 'time_series', 
    #                     '__globals__', 'Nvox', '__header__', '__version__']
    to_check = sio.loadmat(mat_file)
    nvox = to_check['Nvox'][0]
    nb_runs = to_check['time_series'].shape[2] # has shape (time, rois, nb_runs)
    assert nb_runs == 4

    # make a dict for nvox
    check_nvox = {}
    for idx, roi in enumerate(to_check['all_rois']):
        k, _ = osp.splitext(osp.basename(roi[0][0]))
        check_nvox[k] = nvox[idx]

    # make a dict for signals
    arielle_runs = []
    for run in range(nb_runs):
        check_signals = {}
        for idx, roi in enumerate(to_check['all_rois']):
            k  = osp.splitext(osp.basename(roi[0][0]))[0]
            check_signals[k] = to_check['time_series'][:,idx,run]
        arielle_runs.append(check_signals)

    return check_nvox, arielle_runs

def _get_basedir():

    # check on which machine we are working:
    if osp.isdir(BASEDIR_JB):
        base_dir = BASEDIR_JB
    elif osp.isdir(BASEDIR_NX):
        base_dir = BASEDIR_NX

    return base_dir

def _get_smpace_processing_data():

    base_dir = _get_basedir()
    params = get_params(base_dir)
    # modify params to not write signals
    params['analysis']['write_signals'] = False
    params['analysis']['apply_sess_mask'] = False
    params['data']['nb_sess'] = 1
    nb_runs = 4

    # get signals, info 
    subjs_info = process_all(base_dir, params=params, verbose=False)

    jb_runs = []
    for run in range(nb_runs):
        run_str = 'run{:02d}'.format(run+1)
        jb_runs.append(subjs_info['sub01']['sess01'][run_str])
        # has key: ['info', 'arr_sig', 'labels_sig', 'issues', 'arr_sig_f', ...etc]

    return jb_runs

def test_process_all():
    
    check_nvox, arielle_runs = _get_matfile_data()
    jb_runs = _get_smpace_processing_data()

    # check the number of voxels for run 0 (same across runs):
    run0 = 0
    for k in jb_runs[run0]['labels_sig']:
        assert check_nvox[k] == jb_runs[run0]['info'][k]
    for k in jb_runs[run0]['issues']:
        assert check_nvox[k] == 0, "nb vox at {} is {}".format(k, check_nvox[k])

    # check the signals 
    for run in range(len(jb_runs)):
        print('.')
        for idx,k in enumerate(jb_runs[run]['labels_sig']):
            assert_allclose(arielle_runs[run][k], jb_runs[run]['arr_sig'][:,idx])

    #jb_nvox = np.asarray([(k,info_sig[k]) for k in sorted(info_sig.keys())])
    #assert_array_equal(np.asarray(arrielle_nvox), np.asarray(jb_nvox))

def test_csf_extracted_signal():


    check_nvox, arielle_runs = _get_matfile_data()
    

    # print(dir(suf))
    base_dir = _get_basedir()
    params = get_params(base_dir)
    dlayo = params['layout']
    nb_run = params['data']['nb_run']

    sub, sess = 1,1
    runs_dir = osp.join(base_dir, dlayo['dir']['sub+'].format(sub), 
                                  dlayo['dir']['sess+'].format(sess),
                                  dlayo['dir']['runs'])
    csf_dir = osp.join(runs_dir, dlayo['csf']['dir'])
    #print(csf_dir, ' + pat ', dlayo['csf']['roi_mask'])
    csf_file = gb.glob(osp.join(csf_dir, dlayo['csf']['roi_mask']))
    if not csf_file: print("glob empty: {} {}".format(
                                csf_dir, dlayo['csf']['roi_mask']))
    csf_file = suf._check_glob_res(csf_file, ensure=1, files_only=True)

    dir_smooth_imgs = osp.join(runs_dir, dlayo['dir']['smooth'])   

    pat_imgs_files = dlayo['pat']['sub+sess+run+']+"*.nii*"
    runs_pat = [pat_imgs_files.format(sub, sess, run_idx) \
                                        for run_idx in range(1, nb_run+1)]

    runs = [gb.glob(osp.join(dir_smooth_imgs, pat)) for pat in runs_pat]
    for run in runs: run.sort()

    #print("csf_file : ", csf_file)
    #print("csf_dir : ", csf_dir)
    #print("dir_smooth_imgs : ", dir_smooth_imgs)
    #for k in arielle_runs[0].keys():
    #    print(k)

    for idx, run in enumerate(runs):
        run_4d = concat_niimgs(run, ensure_ndim=4)
        #--- get CSF
        csf_arr, csf_labs = ucr.extract_roi_run(csf_dir, csf_file, 
                                run_4d, standardize=False, check_lengh=196, verbose=False)

        mat_arr = arielle_runs[idx]['csf_map_0.93_erode']
        # todo : 'wm_map_0.99_erode'
        assert_allclose(mat_arr, csf_arr)
        #print(csf_labs)
        #print('csf_map_0.93_erode' in arielle_runs[idx].keys())

    pass


