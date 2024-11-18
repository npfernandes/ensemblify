"""Back-end for the Ensemblify Command Line Interface."""

# IMPORTS
## Standard Library Imports
import argparse
import sys

# CONSTANTS
NO_ARGS_ERROR_MSG = '''
Error: Missing required arguments.
Usage: ensemblify {generation,conversion,analysis,reweighting,clash_checking} [module options]
Run \'ensemblify help\' for more information.

'''

ENSEMBLIFY_HELP_MSG = '''
usage: ensemblify {generation,conversion,analysis,reweighting,clash_checking,pipeline,help} [module options]

Command-line tool to access the modules of the Ensemblify Python library.

positional arguments:

    help (h)                   Show this message and exit.
    generation (g, gen)        Access the generation module.
    conversion (c, con)        Access the conversion module.
    analysis (a, ana)          Access the analysis module.
    reweighting (r, rew)       Access the reweighting module.
    clash_checking (cch)       Access the clash checking module.
    pipeline (pip)             Access the pipeline module.

'''

# CLASSES
class CustomHelpFormatter(argparse.HelpFormatter):
    """Helper class derived from argparse.HelpFormatter to make our help menus cleaner.
    Looks hackish, namely the almost full cloning of the _format_action method but it is
    unavoidable given the way argparse is built.
    """
    def __init__(self, prog, max_help_position=4,width=80):
        super().__init__(prog, max_help_position,width)

    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        else:
            return ', '.join(action.option_strings)

    def _format_action(self, action):
        # determine the required width and the entry label
        help_position = min(self._action_max_length + 2,
                            self._max_help_position)
        help_width = max(self._width - help_position, 11)
        action_width = help_position - self._current_indent - 2
        action_header = self._format_action_invocation(action)

        # no help; start on same line and add a final newline
        if not action.help:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup

        # short action name; start on the same line and pad two spaces
        elif len(action_header) <= action_width:
            tup = self._current_indent, '', action_width, action_header
            action_header = '%*s%-*s  ' % tup
            indent_first = 0

        # long action name; start on the next line
        else:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = help_position

        # collect the pieces of the action help
        parts = [action_header]

        # if there was help for the action, add lines of help text
        if action.help and action.help.strip():
            help_text = self._expand_help(action)
            if help_text:
                help_lines = self._split_lines(help_text, help_width)
                parts.append('%*s%s\n' % (indent_first, '', help_lines[0]))
                for line in help_lines[1:]:
                    parts.append('%*s%s\n' % (help_position, '', line))

        # or add a newline if the description doesn't end with one
        elif not action_header.endswith('\n'):
            parts.append('\n')

        # if there are any sub-actions, add their help as well
        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))

        if action_header.endswith(('gen)\n','con)\n','ana)\n','rew)\n','cch)\n','pip)\n')):
            metavar_msg = '    ' + parts[0].strip()
            help_msg = parts[1].lstrip()

            # Extracting the access word index
            access_word_index = help_msg.index('Access')

            # Creating aligned help message
            aligned_help_msg = help_msg[:access_word_index] +\
                               help_msg[access_word_index:].rjust(len(help_msg)\
                                                                   - access_word_index)

            parts = '{:<30} {}'.format(metavar_msg, aligned_help_msg)

        # return a single string
        return self._join_parts(parts)

# FUNCTIONS
def main():
    # Create the initial argument parser to capture the module
    initial_parser = argparse.ArgumentParser(description=('Command-line tool to access various '
                                                          'modules of the Ensemblify Python '
                                                          'library.'),
                                             usage=('ensemblify {generation,conversion,analysis,'
                                                    'reweighting,clash_checking,pipeline}'
                                                    ' [module options]'),
                                             add_help=False) # required

    initial_parser.add_argument('module', choices=['generation', 'g', 'gen',
                                                   'conversion', 'c', 'con',
                                                   'analysis', 'a', 'ana',
                                                   'reweighting', 'r', 'rew',
                                                   'clash_checking', 'cch',
                                                   'pipeline', 'ppl',
                                                   'h', 'help'])

    # Print error message if no arguments are provided
    if len(sys.argv) == 1:
        print(NO_ARGS_ERROR_MSG)
        sys.exit(1)

    # Parse only the module first
    args, remaining_args = initial_parser.parse_known_args()

    # Now create the full parser with subparsers for each module
    parser = argparse.ArgumentParser(description=('Command-line tool to access the various '
                                                  'modules of the Ensemblify Python library.'),
                                     usage='ensemblify',
                                     formatter_class=CustomHelpFormatter)

    # Create subparsers for the modules
    subparsers = parser.add_subparsers(dest='module',
                                       required=True,
                                       metavar='')

    # Expose the help option since initial parser has no help menu
    subparsers.add_parser(name='help',
                          aliases=['h'])

    # Subparser for the 'pipeline' module with aliases
    parser_pipeline = subparsers.add_parser(name='pipeline',
                                            help='Access the pipeline module',
                                            aliases=['ppl'],
                                            usage='ensemblify {pipeline,pip} [options]',
                                            description=('The pipeline module of the '
                                                         'Ensemblify Python library.'),
                                            formatter_class=CustomHelpFormatter)

    parser_pipeline.add_argument('-p', '--parameters',
                                 type=str, required=True, metavar='',
                                 help='Path to parameters file (.yaml).')

    parser_pipeline.add_argument('-a',
                                 action='store_true', default=False,
                                 help=('(Optional) Whether to perform the analysis of the '
                                       'ensemble. Defaults to False.'))

    parser_pipeline.add_argument('-e',
                                 default=None, type=str, metavar='',
                                 help='(Optional) Path to experimental SAXS data file (.dat).')

    # Subparser for the 'clash_checking' module with aliases
    parser_clash_check = subparsers.add_parser(name='clash_checking',
                                               help='Access the clash checking module',
                                               aliases=['cch'],
                                               usage='ensemblify {clash_checking,cch} [options]',
                                               description=('The clash checking module of the '
                                                            'Ensemblify Python library.'),
                                               formatter_class=CustomHelpFormatter)

    parser_clash_check.add_argument('-e', '--ensembledir',
                                    type=str, required=True, metavar='',
                                    help=('Path to directory where ensemble .pdb structures '
                                          'are stored.'))

    parser_clash_check.add_argument('-s',
                                    default=None, type=str, metavar='',
                                    help=('(Optional) Path to file (.yaml) with sampling targets: '
                                          'mapping of chain letters to residue ranges.'))

    parser_clash_check.add_argument('-i',
                                    default=None, type=str, metavar='',
                                    help=('(Optional) Path to input structure (.pdb) used to '
                                          'generate the ensemble.'))

    # Subparser for the 'generation' module with aliases
    parser_generation = subparsers.add_parser(name='generation',
                                              help='Access the generation module',
                                              aliases=['g', 'gen'],
                                              usage='ensemblify {generation,gen,g} [options]',
                                              description=('The generation module of the '
                                                           'Ensemblify Python library.'),
                                              formatter_class=CustomHelpFormatter)

    parser_generation.add_argument('-p', '--parameters',
                                   type=str, required=True, metavar='',
                                   help='Path to parameters file (.yaml).')

    # Subparser for the 'conversion' module with aliases
    parser_conversion = subparsers.add_parser(name='conversion',
                                              help='Access the conversion module',
                                              aliases=['c', 'con'],
                                              usage='ensemblify {conversion,con,c} [options]',
                                              description=('The conversion module of the '
                                                           'Ensemblify Python library.'),
                                              formatter_class=CustomHelpFormatter)

    parser_conversion.add_argument('-j', '--jobname',
                                   required=True, type=str, metavar='',
                                   help='Name for created trajectory file (.xtc).')

    parser_conversion.add_argument('-e', '--ensembledir',
                                   required=True, type=str,  metavar='',
                                   help=('Path to directory where ensemble files (.pdb) are '
                                         'located.'))

    parser_conversion.add_argument('-t', '--trajectorydir',
                                   required=True, type=str,  metavar='',
                                   help=('Path to directory where trajectory file (.xtc) will be '
                                         'created.'))

    parser_conversion.add_argument('-s',
                                   default=10000, type=int,  metavar='',
                                   help=('(Optional) Number of frames of created trajectory file '
                                         '(.xtc).'))

    # Subparser for the 'analysis' module with aliases
    parser_analysis = subparsers.add_parser(name='analysis',
                                            help='Access the analysis module',
                                            aliases=['a', 'ana'],
                                            usage='ensemblify {analysis,ana,a} [options]',
                                            description=('The analysis module of the Ensemblify '
                                                         'Python library.'),
                                            formatter_class=CustomHelpFormatter)

    parser_analysis.add_argument('-trj', '--trajectory',
                                 nargs='+', required=True, type=str,  metavar='',
                                 help='Path(s) to trajectory file(s) (.xtc).')

    parser_analysis.add_argument('-top', '--topology',
                                 nargs='+', required=True, type=str,  metavar='',
                                 help='Path(s) to topology file(s) (.pdb).')

    parser_analysis.add_argument('-tid', '--trajectoryid',
                                 nargs='+', required=True, type=str,  metavar='',
                                 help='Prefix identifier(s) for trajectory file(s).')

    parser_analysis.add_argument('-out',
                                 type=str,  metavar='',
                                 help='(Optional) Path to output directory.')

    parser_analysis.add_argument('-rma',
                                 action ='store_false',
                                 help=('(Optional) Whether to calculate a dihedral angles '
                                       'matrix. Defaults to True.'))

    parser_analysis.add_argument('-dm',
                                 action ='store_false',
                                 help=('(Optional) Whether to calculate a distance matrix. '
                                       'Defaults to True.'))

    parser_analysis.add_argument('-cm',
                                 action ='store_false',
                                 help=('(Optional) Whether to calculate a contact matrix. '
                                       'Defaults to True.'))

    parser_analysis.add_argument('-ssf',
                                 action ='store_false',
                                 help=('(Optional) Whether to calculate a secondary structure '
                                       'assignment frequency matrix. Defaults to True.'))

    parser_analysis.add_argument('-rg',
                                 action ='store_false',
                                 help=('(Optional) Whether to calculate a radius of gyration '
                                       'distribution. Defaults to True.'))
    
    parser_analysis.add_argument('-dmax',
                                 action ='store_false',
                                 help=('(Optional) Whether to calculate a maximum distance '
                                       'distribution. Defaults to True.'))

    parser_analysis.add_argument('-eed',
                                 action ='store_false',
                                 help=('(Optional) Whether to calculate an end-to-end distance '
                                       'distribution. Defaults to True.'))

    parser_analysis.add_argument('-cmdist',
                                 nargs='+', type=str, metavar='',
                                 help=('(Optional) Pair(s) of MDAnalysis selection strings for '
                                       'which to calculate a center of mass distance '
                                       'distribution. Defaults to None.'))

    parser_analysis.add_argument('-colors',
                                 nargs='+', type=str,  metavar='',
                                 help=('(Optional) List of color hexcodes to use, one for each '
                                       'analyzed trajectory.'))

    # Subparser for the 'reweighting' module with aliases
    parser_reweighting = subparsers.add_parser(name='reweighting',
                                               aliases=['r', 'rew'],
                                               help='Access the reweighting module',
                                               usage='ensemblify {reweighting,rew,r} [options]',
                                               description=('The reweighting module of the '
                                                            'Ensemblify Python library.'),
                                               formatter_class=CustomHelpFormatter)

    parser_reweighting.add_argument('-trj','--trajectory',
                                    required=True, type=str,  metavar='',
                                    help='Path to trajectory file (.xtc).')

    parser_reweighting.add_argument('-top', '--topology',
                                    required=True, type=str,  metavar='',
                                    help='Path to topology file (.pdb).')

    parser_reweighting.add_argument('-tid', '--trajectoryid',
                                    required=True, type=str,  metavar='',
                                    help='Prefix identifier for trajectory file.')

    parser_reweighting.add_argument('-exp', '--expdata',
                                    required=True, type=str,
                                    help='Path to experimental SAXS data file (.dat).')

    parser_reweighting.add_argument('-out',
                                    type=str,  metavar='',
                                    help=('(Optional) Path to output directory. Defaults '
                                          'to current working directory.'))

    parser_reweighting.add_argument('-tht',
                                    nargs='+', type=int,  metavar='',
                                    help=('(Optional) List of values to try as the theta '
                                          'parameter in BME.'))

    parser_reweighting.add_argument('-cm',
                                    type=str,
                                    help=('(Optional) Path to calculated contact matrix file '
                                          '(.csv).'))

    parser_reweighting.add_argument('-dm',
                                    type=str,
                                    help=('(Optional) Path to calculated distance matrix file '
                                          '(.csv).'))

    parser_reweighting.add_argument('-ss',
                                    type=str,
                                    help=('(Optional) Path to calculated secondary structure '
                                          'frequency matrix file (.csv).'))

    parser_reweighting.add_argument('-m',
                                    type=str,
                                    help=('(Optional) Path to calculated structural metrics file '
                                          '(.csv).'))

    parser_reweighting.add_argument('-crg',
                                    action ='store_true',
                                    help=('(Optional) Whether to calculate and compare '
                                          'uniform/reweighted radius of gyration distributions. '
                                          'Defaults to False.'))

    parser_reweighting.add_argument('-cdmax',
                                    action ='store_true',
                                    help=('(Optional) Calculate and compare uniform/reweighted '
                                          'maximum distance distributions. Defaults to False.'))

    parser_reweighting.add_argument('-ceed',
                                    action ='store_true',
                                    help=('(Optional) Calculate and compare uniform/reweighted '
                                          'end-to-end distance distributions. Defaults to False.'))

    parser_reweighting.add_argument('-ccmdist',
                                    nargs='+', type=str, metavar='',
                                    help=('(Optional) Pair(s) of MDAnalysis selection strings for '
                                          'which to calculate and compare uniform/reweighted '
                                          'center of mass distance distributions. Defaults to '
                                          'None.'))

    # Now parse the remaining arguments with the full parser
    full_args = parser.parse_args([args.module] + remaining_args)

    # Handle the different modules based on the parsed arguments
    if full_args.module in ['help','h']:
        print(ENSEMBLIFY_HELP_MSG)

    elif full_args.module in ['pipeline','ppl']:
        from ensemblify.pipeline import ensemblify_pipeline
        ensemblify_pipeline(parameters=full_args.parameters,
                            analysis=full_args.a,
                            exp_data=full_args.e,)

    elif full_args.module in ['clash_checking', 'cch']:
        from ensemblify.clash_checking import check_steric_clashes
        check_steric_clashes(ensemble_dir=full_args.ensembledir,
                             sampling_targets=full_args.s,
                             input_structure=full_args.i,)

    elif full_args.module in ['generation', 'g', 'gen']:
        from ensemblify.generation import generate_ensemble
        generate_ensemble(parameters_path=full_args.parameters,)

    elif full_args.module in ['conversion', 'c', 'con']:
        from ensemblify.conversion import ensemble2traj
        ensemble2traj(job_name=full_args.jobname,
                      ensemble_dir=full_args.ensembledir,
                      trajectory_dir=full_args.trajectorydir,
                      trajectory_size=full_args.s,)

    elif full_args.module in ['analysis', 'a', 'ana']:
        from ensemblify.analysis import analyze_trajectory
        analyze_trajectory(trajectories=full_args.trajectory,
                           topologies=full_args.topology,
                           trajectory_ids=full_args.trajectoryid,
                           output_directory=full_args.out,
                           ramachandran_data=full_args.rma,
                           distancematrices=full_args.dm,
                           contactmatrices=full_args.cm,
                           ssfrequencies=full_args.ssf,
                           rg=full_args.rg,
                           dmax=full_args.dmax,
                           eed=full_args.eed,
                           cm_dist=full_args.cmdist,
                           color_palette=full_args.colors,)

    elif full_args.module in ['reweighting', 'r', 'rew']:
        from ensemblify.reweighting import reweight_ensemble
        reweight_ensemble(trajectory=full_args.trajectory,
                          topology=full_args.topology,
                          trajectory_id=full_args.trajectoryid,
                          exp_saxs_data=full_args.expdata,
                          output_dir=full_args.out,
                          thetas=full_args.tht,
                          calculated_cmatrix=full_args.cm,
                          calculated_dmatrix=full_args.dm,
                          calculated_ss_frequency=full_args.ss,
                          calculated_metrics_data=full_args.m,
                          compare_rg=full_args.crg,
                          compare_dmax=full_args.cdmax,
                          compare_eed=full_args.ceed,
                          compare_cmdist=full_args.ccmdist,)
