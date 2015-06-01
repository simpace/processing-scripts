from __future__ import print_function, division
from nose.tools import raises

import numpy as np
from numpy.testing import assert_allclose
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

    # get signals from info
    subjs_info = process_all(base_dir, params=params, verbose=False)
    signals = subjs_info['sub01']['sess01']['run01']['signals'] 
    sizes = signals['info']
    signals = signals['arr_sig']
    labels = signals['labels_sig']


    # compare to Arielle's
    mat_file = osp.join(osp.realpath(__file__), "sub01_session_1_raw_ROI_timeseries.mat")
    to_check = np.load_mat(mat_file)

    # to_check.keys() == ['all_rois', 'time_series', 
    #                     '__globals__', 'Nvox', '__header__', '__version__']
    nvox = to_check['Nvox'][0]
    for roi in to_check['all_rois']:
        k, _ = osp.splitext(osp.basename(roi[0][0]))
        if k in 


    


