"""Unit tests for the ensemblify.generation module."""

# IMPORTS
## Standard Library Imports
#import pytest

## Local Imports
from ensemblify.generation import generate_ensemble
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

def test_generate_ensemble():
    input_params = 'examples/input_parameters/generate_Hst5.yaml'
    valid_pdbs_dir = generate_ensemble(input_params)
    assert valid_pdbs_dir == './examples/Hst5_test/ensemble/valid_pdbs', valid_pdbs_dir

if __name__ == '__main__':
    test_generate_ensemble()
