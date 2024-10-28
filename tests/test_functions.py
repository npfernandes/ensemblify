"""Unit tests for the Ensemblify Python Library."""

# IMPORTS
## Standard Library Imports
#import pytest

## Local Imports
from ensemblify.generation.ensemble_utils.functions import setup_pose, derive_constraint_targets

def test_derive_constraint_targets():
    pose = setup_pose('AAAAAAAAAAAAAAAAAAAAAAAA') # 24 residues

    result_1 = derive_constraint_targets(pose,{'A' : (('MC',(6,8),'all','TRIPEPTIDE'),
                                                      ('MC',(12,16),'all','TRIPEPTIDE'),
                                                      ('MC',(20,22),'all','TRIPEPTIDE'))})
    assert result_1 == ((1,5),(9,11),(17,19),(23,24)), result_1

    result_2 = derive_constraint_targets(pose,{'A' : (('MC',(2,4),'all','TRIPEPTIDE'),
                                                      ('MC',(6,8),'all','TRIPEPTIDE'))})
    assert result_2 == ((1,1),(5,5),(9,24)), result_2

    result_3 = derive_constraint_targets(pose,{'A' : (('MC',(1,3),'all','TRIPEPTIDE'),
                                                      ('MC',(7,24),'all','TRIPEPTIDE'))})
    assert result_3 == ((4,6),), result_3

    result_4 = derive_constraint_targets(pose,{'A' : (('MC',(1,4),'all','TRIPEPTIDE'),)})
    assert result_4 == ((5,24),), result_4

    result_5 = derive_constraint_targets(pose,{'A' : (('MC',(7,24),'all','TRIPEPTIDE'),)})
    assert result_5 == ((1,6),), result_5

    result_6 = derive_constraint_targets(pose,{'A' : (('MC',(2,8),'all','TRIPEPTIDE'),)})
    assert result_6 == ((1,1),(9,24)), result_6

    result_7 = derive_constraint_targets(pose,{'A' : (('MC',(1,2),'all','TRIPEPTIDE'),
                                                      ('MC',(4,6),'all','TRIPEPTIDE'),
                                                      ('MC',(8,24),'all','TRIPEPTIDE'))})
    assert result_7 == ((3,3),(7,7)), result_7

    result_8 = derive_constraint_targets(pose,{'A' : (('MC',(1,4),'all','TRIPEPTIDE'),
                                                      ('MC',(8,10),'all','TRIPEPTIDE'))})
    assert result_8 == ((5,7),(11,24)), result_8

    result_9 = derive_constraint_targets(pose,{'A' : (('MC',(6,8),'all','TRIPEPTIDE'),
                                                      ('MC',(10,24),'all','TRIPEPTIDE'))})
    assert result_9 == ((1,5),(9,9)), result_9

# Clash checking

# if __name__ == '__main__':
#     from ensemblify import update_config
#     update_config({'FASPR_PATH':'/home/tiagogomes/software/FASPR-master/FASPR',
#                    'PULCHRA_PATH':'/home/tiagogomes/software/pulchra-master/pulchra_CHANGED'})
    
    
#     SAMPLING_TARGETS = {'A' : [ [ 'MC', [1,24], 'coil', 'TRIPEPTIDE' ]]}
#     INPUT_STRUCTURE = '/home/tiagogomes/Desktop/projects/nuno_fernandes/proteins_plus_saxs/starting_structures_atomistic/IDPs/Hst5.pdb'
#     ENSEMBLE_DIR = '/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/ENSEMBLES/Hst5/ENSEMBLE_Rechecked'

#     check_steric_clashes(ensemble_dir='/home/tiagogomes/Downloads/protein_pool',
#                          sampling_targets=None,
#                          input_structure=None)

# Analyze trajectory

# if __name__ == '__main__':
#     analyze_trajectory(['/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/Hst5/Hst5_trajectory.xtc'],
#                        ['/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/Hst5/Hst5_top.pdb'],
#                        ['Hst5'],
#                        output_directory='/home/tiagogomes/Desktop/projects/nuno_fernandes/NProtein_sarscov2/NProtein_365_TetramerClosed_Ensemble/testing_hst5_anal',
#                        ramachandran_data=False,
#                        contactmatrices=False,
#                        distancematrices=False,
#                        ssfrequencies=False,
#                     #    rg=False,
#                     #    dmax=False,
#                     #    eed=False
#                     )

#     # analyze_trajectory(['/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_from_AlphaFold/TRAJECTORIES/LysP7951/LysP7951_trajectory.xtc'],
#     #                    ['/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_from_AlphaFold/TRAJECTORIES/LysP7951/LysP7951_top.pdb'],
#     #                    ['LysP7951'],
#     #                    output_directory='/home/tiagogomes/Desktop/projects/nuno_fernandes/NProtein_sarscov2/NProtein_365_TetramerClosed_Ensemble/testing_lysp7951_anal',
#     #                    ramachandran_data=False,
#     #                     contactmatrices=False,
#     #                    #distancematrices=False,
#     #                    ssassignments=False,
#     #                    rg=False,
#     #                    dmax=False,
#     #                    eed=False)