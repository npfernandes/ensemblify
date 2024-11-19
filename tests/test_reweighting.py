# """Unit tests for the ensemblify.reweighting module."""

# # IMPORTS
# ## Standard Library Imports
# #import pytest

# ## Local Imports


# if __name__ == '__main__':
#     from ensemblify import update_config

#     update_config({'PEPSI_SAXS_PATH': '/home/tiagogomes/software/Pepsi-SAXS',
#                    'BIFT_PATH': '/home/tiagogomes/software/bift'})

#     reweight_ensemble(trajectory='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/Hst5/Hst5_trajectory.xtc',
#                       topology='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/Hst5/Hst5_top.pdb',
#                       trajectory_id='Hst5',
#                       exp_saxs_data='/home/tiagogomes/Desktop/projects/nuno_fernandes/proteins_plus_saxs/SAXS/bift_Hst5.dat',
#                       output_dir='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/REWEIGHTING',
#                       calculated_cmatrix='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/Hst5/Hst5_contact_matrix.csv',
#                       calculated_dmatrix='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/Hst5/Hst5_distance_matrix.csv',
#                       calculated_ss_frequency='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/Hst5/Hst5_ss_frequency.csv',
#                       calculated_metrics_data='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/Hst5/Hst5_structural_metrics.csv'
#                     )

#     # reweight_ensemble(trajectory='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/THB_C2/THB_C2_trajectory.xtc',
#     #                   topology='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/THB_C2/THB_C2_top.pdb',
#     #                   trajectory_id='THB_C2',
#     #                   exp_saxs_data='/home/tiagogomes/Desktop/projects/nuno_fernandes/proteins_plus_saxs/SAXS/bift_THB_C2.dat',
#     #                   output_dir='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/REWEIGHTING',
#     #                   calculated_cmatrix='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/THB_C2/THB_C2_contact_matrix.csv',
#     #                   calculated_dmatrix='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/THB_C2/THB_C2_distance_matrix.csv',
#     #                   calculated_ss_frequency='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/THB_C2/THB_C2_ss_frequency.csv',
#     #                   calculated_metrics_data='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORY_ANALYSIS/THB_C2/THB_C2_structural_metrics.csv'
#     #                 )
