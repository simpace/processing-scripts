from __future__ import print_function, division
from nose.tools import raises

from ..smpce_data_to_corr import process_all

@raises(IOError)
def test_process_all_assert():
    process_all('1')


