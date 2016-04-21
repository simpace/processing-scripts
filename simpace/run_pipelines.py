
import os
import json

base_dir = '/home/despo/simpace/rename_files/derivatives'
subjects = ['1', '2']
traj = 'trajA'
nsess = [13, 6]
params = range(0, 5)
print params

for idx, iss in enumerate(subjects):

    sessions = str(range(1, nsess[idx]+1))
    sessions = sessions.replace(' ', '')

    nuisance_file = 'analysis_parameters_extract_LFF.json'

    # for iparam in params:
    #     with open(nuisance_file) as fnuis:
    #         dnuis = json.load(fnuis)
    #         dnuis['data_prefix']['dir'] = 'NR_' + str(iparam)
    #         json.dump(fnuis, dnuis)

    command = 'python smpce_data_to_corr_mod.py {0} {1} {2} {3} {4}'.format(base_dir, iss, traj, sessions, nuisance_file)
    os.system(command)

    # for iparam in params:
    #
    #     nuisance_file = 'nuisance_parameters_' + str(iparam) + '.json'
    #
    #     command = 'python smpce_data_to_corr_mod.py {0} {1} {2} {3} {4}'.format(base_dir, iss, traj, sessions, nuisance_file)
    #     os.system(command)

