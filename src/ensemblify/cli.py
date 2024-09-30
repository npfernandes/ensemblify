"""Back-end for the Ensemblify Command Line Interface."""

# IMPORTS
## Standard Library Imports
import argparse
import sys

## Third Party Imports
from ensemblify import generate_ensemble,ensemble2traj,analyze_trajectory,reweigh_ensemble

# CLASSES
class CustomHelpFormatter(argparse.HelpFormatter):
    """Helper class derived from argparse.HelpFormatter to make our help menus cleaner."""
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

        if action_header.endswith(('gen)\n','con)\n','ana)\n','rew)\n')):
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
def read_dict():
    pass

def read_list():
    pass


def main():
    # Create top-level argument parser
    parser = argparse.ArgumentParser(description='Command-line tool to access the various modules of the Ensemblify Python library.',
                                     usage='ensemblify {generation,conversion,analysis,reweighting} [module options]',
                                     formatter_class=CustomHelpFormatter)
    
    # Create subparsers for the modules
    subparsers = parser.add_subparsers(dest='module',
                                       required=True,
                                       metavar='')

    # Subparser for the 'generation' module with aliases
    parser_generation = subparsers.add_parser(name='generation',
                                              help='Access the generation module',
                                              aliases=['g', 'gen'],
                                              usage='ensemblify {generation,gen,g} [options]',
                                              description='The generation module of the Ensemblify Python library.',
                                              formatter_class=CustomHelpFormatter)
    
    parser_generation.add_argument('-p','--parameters', type=str, required=True, help='Path to parameters file in .yaml format.', metavar='PARAMS')

    # Subparser for the 'conversion' module with aliases
    parser_conversion = subparsers.add_parser(name='conversion',
                                              help='Access the conversion module',
                                              aliases=['c', 'con'],
                                              usage='ensemblify {conversion,con,c} [options]',
                                              description='The conversion module of the Ensemblify Python library.',
                                              formatter_class=CustomHelpFormatter)

    parser_conversion.add_argument('-n','--jobname', required=True, type=str, help='Name for created trajectory file (.xtc).', metavar='')
    parser_conversion.add_argument('-i','--ensembledir', required=True, type=str, help='Path to directory where ensemble files (.pdb) are located.', metavar='')
    parser_conversion.add_argument('-o','--trajectorydir', required=True, type=str, help='Path to directory where trajectory file (.xtc) will be created.', metavar='')
    parser_conversion.add_argument('-s','--size', type=int, help='(Optional) Number of frames of created trajectory file (.xtc).', metavar='')

    # Subparser for the 'analysis' module with aliases
    parser_analysis = subparsers.add_parser(name='analysis',
                                            help='Access the analysis module',
                                            aliases=['a', 'ana'],
                                            usage='ensemblify {analysis,ana,a} [options]',
                                            description='The analysis module of the Ensemblify Python library.',
                                            formatter_class=CustomHelpFormatter)
    
    parser_analysis.add_argument('-tj','--trajectory', required=True, type=list[str], help='Path(s) to trajectory file(s) (.xtc).', metavar='')
    parser_analysis.add_argument('-tp','--topology', required=True, type=list[str], help='Path(s) to topology file(s) (.pdb).', metavar='')
    parser_analysis.add_argument('-ti','--trajectoryid', required=True, type=list[str], help='Prefix identifier(s) for trajectory file(s).', metavar='')

    parser_analysis.add_argument('--outputdir', type=str, help='(Optional) Path to output directory.', metavar='')
    parser_analysis.add_argument('--ramadata', action ='store_true', help='(Optional) Calculate a dihedral angles matrix.')
    parser_analysis.add_argument('--distancematrix', action ='store_true', help='(Optional) Calculate a distance matrix.')
    parser_analysis.add_argument('--contactmap', action ='store_true', help='(Optional) Calculate a contact map.')
    parser_analysis.add_argument('--ssassign', action ='store_true', help='(Optional) Calculate a secondary structure assignment matrix.')
    parser_analysis.add_argument('--rg', action ='store_true', help='(Optional) Calculate a radius of gyration distribution.')
    parser_analysis.add_argument('--dmax', action ='store_true', help='(Optional) Calculate a maximum distance distribution.')
    parser_analysis.add_argument('--eed', action ='store_true', help='(Optional) Calculate an end-to-end distance distribution.')
    parser_analysis.add_argument('--cmdist', type=dict[str,(str,str)], help='(Optional) Calculate a center of mass distance distribution for each pair of given MDAnalysis selection strings.', metavar='')
    parser_analysis.add_argument('--colors', type=list[str], help='(Optional) List of color hexcodes to use, one for each analyzed trajectory.', metavar='')

    # Subparser for the 'reweighting' module with aliases
    parser_reweighting = subparsers.add_parser(name='reweighting',
                                               aliases=['r', 'rew'],
                                               help='Access the reweighting module',
                                               usage='ensemblify {reweighting,rew,r} [options]',
                                               description='The reweighting module of the Ensemblify Python library.',
                                               formatter_class=CustomHelpFormatter)

    parser_reweighting.add_argument('--trajectory', required=True, type=str, help='Path to trajectory file (.xtc).', metavar='')
    parser_reweighting.add_argument('--topology', required=True, type=str, help='Path to topology file (.pdb).', metavar='')
    parser_reweighting.add_argument('--trajectoryid', required=True, type=str, help='Prefix identifier for trajectory file.', metavar='')
    parser_reweighting.add_argument('--expdata', required=True, type=str, help='Path to experimental SAXS data file (.dat).', metavar='')

    parser_reweighting.add_argument('--outputdir', type=str, help='(Optional) Path to output directory.', metavar='')
    parser_reweighting.add_argument('--thetas', type=list[int], help='(Optional) List of values to try as the theta parameter in BME.', metavar='')
    parser_reweighting.add_argument('--rg', action ='store_true', help='(Optional) Calculate and compare uniform/reweighted radius of gyration distributions.')
    parser_reweighting.add_argument('--dmax', action ='store_true', help='(Optional) Calculate and compare uniform/reweighted maximum distance distributions.')
    parser_reweighting.add_argument('--eed', action ='store_true', help='(Optional) Calculate and compare uniform/reweighted end-to-end distance distributions.')
    parser_reweighting.add_argument('--cmdist', type=dict[str,(str,str)], help='(Optional) Calculate and compare uniform/reweighted center of mass distance distributions for each pair of given MDAnalysis selection strings.', metavar='')

    # Print error message if no arguments are provided
    if len(sys.argv) == 1:
        print('Error: Missing required arguments.\n')
        print('Usage:')
        print('  ensemblify {generation, conversion, analysis, reweighting} [module options]\n')
        print('Run \'ensemblify --help\'  for more information.\n')
        sys.exit(1)

    # Parse the command-line arguments
    args = parser.parse_args()

    # Handle the different modules
    if args.module in ['generation', 'g', 'gen']:
        generate_ensemble(parameters_path=args.PARAMS)

    elif args.module in ['conversion', 'c', 'con']:
        ensemble2traj(job_name=args.jobname,
                      ensemble_dir=args.ensembledir,
                      trajectory_dir=args.trajectorydir)

    elif args.module in ['analysis', 'a', 'ana']:
        analyze_trajectory(trajectories=args.trajectory,
                           topologies=args.topology,
                           trajectory_ids=args.trajectoryid,
                           output_directory=args.outputdir,
                           ramachandran_data=args.ramadata,
                           distancematrices=args.distancematrix,
                           contactmatrices=args.contactmap,
                           ssassignments=args.ssassign,
                           rg=args.rg,
                           dmax=args.dmax,
                           eed=args.eed,
                           cm_dist=args.cmdist,
                           color_palette=args.colors)

    elif args.module in ['reweighting', 'r', 'rew']:
        reweigh_ensemble(trajectory=args.trajectory,
                         topology=args.topology,
                         trajectory_id=args.trajectoryid,
                         exp_saxs_data=args.expdata,
                         output_dir=args.outputdir,
                         thetas=args.thetas,
                         compare_rg=args.rg,
                         compare_dmax=args.dmax,
                         compare_eed=args.eed,
                         compare_cmdist=args.cmdist)
