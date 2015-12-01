from __future__ import print_function

import os.path as osp
import json
import argparse
import smpce_data_to_corr as stc

def main(dbase, json_files, run_default=False, verbose=False):

    """
    """
    # read common params
    default_params = stc.get_params(dbase, verbose=verbose)
    #cwd = osp.dirname(osp.realpath(__file__))
    cwd = dbase 
    if run_default:
        info = stc.process_all("dummy_dbase", default_params, verbose=verbose)
        ftmp = osp.join(cwd,'tmp_' + osp.join('default','log'))
        with open(ftmp, 'w') as this_file: 
            print(info, file=this_file)

    # add what is in the additional json files
    for json_file in json_files: 
        # initialize params to default params
        params = default_params
        # add what is in the additional jason file
        if verbose: print("adding json file {}".format(json_file))

        with open(json_file) as fjson:
            djson = json.load(fjson)
        for djson_k in djson:
            params['layout'][djson_k] = djson[djson_k]
            if verbose: print(djson_k, djson[djson_k])

        info = stc.process_all("dummy_dbase", params, verbose=verbose)

        ftmp = osp.join(cwd,'tmp_' + osp.join(osp.basename(json_file),'log'))

        with open(ftmp, 'w') as this_file: 
            print(info, file=this_file)


#------------------------------------------------------------------------------
# Running from command line
#------------------------------------------------------------------------------

if __name__ == "__main__":
    #
    parser = argparse.ArgumentParser()

    # Positional required arguments
    help_base_dir = "The directory containing sub01/sess??/... and json files: \n" + \
                    "analysis_parameters.json  data_parameters.json  directory_layout.json"
    parser.add_argument("base_directory", help=help_base_dir)

    # Optional arguments

    help_add_json = "an additional json layout file overwritting the one found in the base dir"

    parser.add_argument('-a','--add_jsons', nargs='+', help=help_add_json, required=False)

    parser.add_argument("--run_default", help="run default analsysis", action="store_true")
    
    parser.add_argument("--verbose", help="increase output verbosity", action="store_true")

    args = parser.parse_args()
    base_dir = args.base_directory
    add_jsons = args.add_jsons
    run_default = args.run_default
    verbose = args.verbose
    
    print("launching analyses on :", base_dir)
    assert osp.isdir(base_dir), '{} not a directory'.format(base_dir)

    info = main(base_dir, add_jsons, run_default=run_default, verbose=verbose)
     
    if verbose:
        print("\n------------------- Debug info ------------------ \n")
        print(info)
        print("\n------------------- Debug info ------------------ \n")


