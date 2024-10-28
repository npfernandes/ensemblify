"""View and update Ensemblify's global configuration."""

# IMPORTS
## Standard Library Imports
import os

# CONSTANTS
# Global configuration dictionary
GLOBAL_CONFIG = {
    'USED_DATABASE_COLNAMES': {'OMG1': 'OMG1', # values == column names required in input dbs
                               'OMG2': 'OMG2', # case insensitive, becomes capital letter
                               'OMG3': 'OMG3',
                               'PHI1': 'PHI1',
                               'PHI2': 'PHI2',
                               'PHI3': 'PHI3',
                               'PSI1': 'PSI1',
                               'PSI2': 'PSI2',
                               'PSI3': 'PSI3',
                               'FRAG': 'FRAG'},
    'ALPHA_HELIX_CANON': (-57,-47), # (Phi,Psi); Lehninger principles of biochemistry (2021)
    'BETA_STRAND_CANON': (-135,135), # (Phi,Psi); Wikipedia Beta strand # FIXME
    'FASPR_PATH': os.environ.get('FASPR_PATH'), # get faspr environment variable
    'PULCHRA_PATH': os.environ.get('PULCHRA_PATH'), # get pulchra environment variable
    'PEPSI_SAXS_PATH': os.environ.get('PEPSI_SAXS_PATH'), # get pepsi-saxs environment variable
    'BIFT_PATH': os.environ.get('BIFT_PATH'), # get bift environment variable (optional)
    'PLOTLY_DISPLAY_CONFIG': {'displaylogo': False,
                              'toImageButtonOptions': {'format': 'svg', # defaults to svg download
                                                       'height': None,
                                                       'width': None,
                                                       'scale': 1}} # defaults to bift alias
    }

# FUNCTIONS
def show_config() -> dict:
    """
    Show the current configuration dictionary.

    Returns:
        The Ensemblify global configuration dictionary.
    
    """
    return GLOBAL_CONFIG


def update_config(new_config: dict):
    """
    Update the configuration dictionary with new values.
    
    Any keys in the new configuration that exist in the
    current config are updated.

    Args:
        new_config:
            dictionary with configuration parameters to update.

    """
    GLOBAL_CONFIG.update(new_config)
