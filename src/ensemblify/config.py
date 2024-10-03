"""Auxilliary functions to view and update Ensemblify's global configuration."""

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
    'FASPR_PATH': 'faspr', # defaults to faspr alias
    'PULCHRA_PATH': 'pulchra', # defaults to pulchra alias
    'PEPSI_SAXS_PATH': 'pepsisaxs', # defaults to pepsisaxs alias
    'BIFT_PATH': 'bift', # defaults to bift alias
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
    
    Example:
        >>> import ensemblify as ey
        >>> ey.show_config()
        {'USED_DATABASE_COLNAMES': {'OM1': 'OMG1',
                                    'OM2': 'OMG2',
                                    'OM3': 'OMG3',
                                    'PH1': 'PHI1',
                                    'PH2': 'PHI2',
                                    'PH3': 'PHI3',
                                    'PS1': 'PSI1',
                                    'PS2': 'PSI2',
                                    'PS3': 'PSI3',
                                    'FRG': 'FRAG'}
         'ALPHA_HELIX_CANON': (-57,-47),
         'BETA_STRAND_CANON': (-135,135),
         'FASPR_PATH': 'faspr',
         'PULCHRA_PATH': 'pulchra',
         'PEPSI_SAXS_PATH': 'pepsi-saxs',
         'BIFT_PATH': 'bift',
         'PLOTLY_DISPLAY_CONFIG': {'displaylogo': False,
                                   'toImageButtonOptions': {'format': 'svg',
                                   'height': None,
                                   'width': None,
                                   'scale': 1}}}
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

    Example:
        >>> import ensemblify as ey
        >>> ey.show_config()
        {... , 'FASPR_PATH': 'faspr', ...}
        >>> ey.update_config({'FASPR_PATH': 'your_faspr_path_here'})
        >>> ey.show_config()
        {... , 'FASPR_PATH': 'your_faspr_path_here', ...}
    """
    GLOBAL_CONFIG.update(new_config)
