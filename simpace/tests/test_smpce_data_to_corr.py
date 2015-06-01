from __future__ import print_function, division
from nose.tools import raises

import numpy as np
import scipy.io as sio
from numpy.testing import assert_allclose #, assert_array_equal
import os.path as osp
from ..smpce_data_to_corr import get_params, process_all

BASEDIR_JB = '/home/jb/data/simpace/data/rename_files'
BASEDIR_NX = '/home/despo/simpace/subject_1_data/rename_files'


@raises(IOError)
def test_process_all_assert():
    process_all('1')

def test_process_all():
    
    # check on which machine we are working:
    if osp.isdir(BASEDIR_JB):
        base_dir = BASEDIR_JB
    elif osp.isdir(BASEDIR_NX):
        base_dir = BASEDIR_NX
    
    params = get_params(base_dir)
    # modify params to not write signals
    params['analysis']['write_signals'] = False
    params['analysis']['apply_sess_mask'] = False
    nb_runs = 4

    # get signals, info 
    subjs_info = process_all(base_dir, params=params, verbose=False)
    #labs_sig = subjs_info['sub01']['sess01']['run01']['signals']['labels_sig']

    jb_runs = []
    for run in range(nb_runs):
        run_str = 'run{:02d}'.format(run+1)
        jb_runs.append(subjs_info['sub01']['sess01'][run_str]['signals'])
        # has key: ['info', 'arr_sig', 'labels_sig', 'issues']

    # compare to Arielle's in the .mat file
    mat_file = osp.join(osp.dirname(osp.realpath(__file__)), "sub01_session_1_raw_ROI_timeseries.mat")
    # to_check.keys() == ['all_rois', 'time_series', 
    #                     '__globals__', 'Nvox', '__header__', '__version__']
    to_check = sio.loadmat(mat_file)
    nvox = to_check['Nvox'][0]
    #time_series = to_check['time_series']

    # make a dict for nvox
    check_nvox = {}
    for idx, roi in enumerate(to_check['all_rois']):
        k, _ = osp.splitext(osp.basename(roi[0][0]))
        check_nvox[k] = nvox[idx]

    # make a dict for signals
    arielle_runs = []
    for run in range(nb_runs):
        check_sign = {}
        for idx, roi in enumerate(to_check['all_rois']):
            k, _ = osp.splitext(osp.basename(roi[0][0]))
            check_sign[k] = to_check['time_series'][:,idx,run]
        arielle_runs.append(check_sign)

    # check the number of voxels:
    for k in jb_runs[0]['labels_sig']:
        assert check_nvox[k] == jb_runs[0]['info'][k]
    for k in jb_runs[0]['issues']:
        assert check_nvox[k] == 0, "nb vox at {} is {}".format(k, check_nvox[k])

    # check the signals 
    for run in range(nb_runs):
        print('.')
        for idx,k in enumerate(jb_runs[run]['labels_sig']):
            assert_allclose(arielle_runs[run][k], jb_runs[run]['arr_sig'][:,idx])

    #jb_nvox = np.asarray([(k,info_sig[k]) for k in sorted(info_sig.keys())])
    #assert_array_equal(np.asarray(arrielle_nvox), np.asarray(jb_nvox))




