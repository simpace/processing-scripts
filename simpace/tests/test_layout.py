from __future__ import division, print_function, absolute_import

# from nose.plugins.attrib import attr

#import json
import layout as lo 
import os.path as osp 
#import smpce_data_to_corr as stc

JSONTEST = '/home/jb/code/simpace/simpace/jsons/essai.json'
DATABASEDIR = '/home/jb/data/simpace/rename_files'

dlayo = lo.get_layout(JSONTEST, dbase = DATABASEDIR)

def test_get_layout():

    assert 'subjects' in dlayo.keys()
    assert 'sessions' in dlayo.keys()
    assert 'runs' in dlayo.keys()
    assert 'run' not in dlayo.keys()

def test_get_key_dict():

    key_dict_expected = {'runs': u'run', 'subjects': u'sub', 'sessions': u'sess'}
    dstate = lo._get_key_dict(dlayo, entities = ["subjects", "sessions", "runs"])
    assert dstate == key_dict_expected

def test_get_pth_globFalse():
    #  get pth with glob = False

    pth, pth_dict = lo._get_pth(dlayo, 'subjects', glob=False)
    assert pth == DATABASEDIR
    assert pth_dict == {}

    pth, pth_dict = lo._get_pth(dlayo, 'sessions', glob=False)
    print("\n",pth)
    expected_pth =  osp.join(DATABASEDIR, 'sub{sub:02d}')
    print('expected_pth: ',expected_pth)
    assert pth == expected_pth
    assert pth_dict == {u'sub': None}

    pth, pth_dict = lo._get_pth(dlayo, 'runs', glob=False)
    print("\n",pth)
    expected_pth = osp.join(DATABASEDIR, 'sub{sub:02d}', 'sess{sess:02d}', 'preproc') 
    print('expected_pth: ',expected_pth)
    assert pth == expected_pth
    assert pth_dict == {u'sess': None, u'sub': None}


def test_get_pth_globTrue():
    #  get pth with glob = True

    pth = lo._get_pth(dlayo, 'subjects', glob=True)
    assert pth == DATABASEDIR

    pth = lo._get_pth(dlayo, 'sessions', glob=True)
    print("\n",pth)
    expected_pth =  osp.join(DATABASEDIR, 'sub*')
    print('expected_pth: ',expected_pth)
    assert pth == expected_pth

    pth = lo._get_pth(dlayo, 'runs', glob=True)
    print("\n",pth)
    expected_pth = osp.join(DATABASEDIR, 'sub*', 'sess*', 'preproc') 
    print('expected_pth: ',expected_pth)
    assert pth == expected_pth



def test_get_pthglb():

    pdir, ddir = lo._get_pthglb(dlayo, 'sessions')



