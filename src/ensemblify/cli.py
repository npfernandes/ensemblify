"""Back-end for the Ensemblify Command Line Interface."""

# IMPORTS
## Standard Library Imports
import argparse
import sys

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
# def read_dict():
#     pass

# def read_list():
#     pass

def main():
    # Create the initial argument parser to capture the module
    initial_parser = argparse.ArgumentParser(description='Command-line tool to access various modules of the Ensemblify Python library.',
                                             usage='ensemblify {generation,conversion,analysis,reweighting} [module options]',
                                             add_help=False) # required

    initial_parser.add_argument('module', choices=['generation', 'g', 'gen', 'conversion', 'c', 'con',
                                                   'analysis', 'a', 'ana', 'reweighting', 'r', 'rew',
                                                   'clash_checking','cch','h','help'])
    
    # Print error message if no arguments are provided
    if len(sys.argv) == 1:
        print('Error: Missing required arguments.\n')
        print('Usage:')
        print('  ensemblify {generation, conversion, analysis, reweighting, clash_checking} [module options]\n')
        print('Run \'ensemblify help\'  for more information.\n')
        sys.exit(1)

    # Parse only the module first
    args, remaining_args = initial_parser.parse_known_args()

    # Now create the full parser with subparsers for each module
    parser = argparse.ArgumentParser(description='Command-line tool to access the various modules of the Ensemblify Python library.',
                                     usage='ensemblify {generation,conversion,analysis,reweighting,clash_checking} [module options]',
                                     formatter_class=CustomHelpFormatter)

    # Create subparsers for the modules
    subparsers = parser.add_subparsers(dest='module',required=True,  metavar='')

    # Expose the help option since initial parser has no help menu
    subparsers.add_parser(name='help',aliases=['h'])

    # Subparser for the 'clash_checking' module with aliases
    parser_clash_check = subparsers.add_parser(name='clash_checking',
                                               help='Access the clash checking module',
                                               aliases=['cch'],
                                               usage='ensemblify {clash_checking,cch} [options]',
                                               description='The clash checking module of the Ensemblify Python library.',
                                               formatter_class=CustomHelpFormatter)

    parser_clash_check.add_argument('-e','--ensemble_dir', type=str, required=True, help='Path to directory where ensemble .pdb structures are stored.',metavar='')
    parser_clash_check.add_argument('-s','--sampling_targets', default=None, type=str, help='(Optional) Path to file (.yaml) with mapping of chains to sampled regions.',metavar='')
    parser_clash_check.add_argument('-i','--input_structure', default=None, type=str, help='(Optional) Path to input structure (.pdb) used to generate the ensemble.',metavar='')

    # Subparser for the 'generation' module with aliases
    parser_generation = subparsers.add_parser(name='generation',
                                              help='Access the generation module',
                                              aliases=['g', 'gen'],
                                              usage='ensemblify {generation,gen,g} [options]',
                                              description='The generation module of the Ensemblify Python library.',
                                              formatter_class=CustomHelpFormatter)

    parser_generation.add_argument('-p','--parameters', type=str, required=True, help='Path to parameters file (.yaml).', metavar='')

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

    parser_conversion.add_argument('-s','--size', default=10000, type=int, help='(Optional) Number of frames of created trajectory file (.xtc).', metavar='')

    # Subparser for the 'analysis' module with aliases
    parser_analysis = subparsers.add_parser(name='analysis',
                                            help='Access the analysis module',
                                            aliases=['a', 'ana'],
                                            usage='ensemblify {analysis,ana,a} [options]',
                                            description='The analysis module of the Ensemblify Python library.',
                                            formatter_class=CustomHelpFormatter)

    parser_analysis.add_argument('--trajectory', required=True, type=list[str], help='Path(s) to trajectory file(s) (.xtc).', metavar='')
    parser_analysis.add_argument('--topology', required=True, type=list[str], help='Path(s) to topology file(s) (.pdb).', metavar='')
    parser_analysis.add_argument('--trajectoryid', required=True, type=list[str], help='Prefix identifier(s) for trajectory file(s).', metavar='')

    parser_analysis.add_argument('--outputdir', type=str, help='(Optional) Path to output directory.', metavar='')
    # parser_analysis.add_argument('--ramadata', action ='store_true', help='(Optional) Calculate a dihedral angles matrix.')
    # parser_analysis.add_argument('--distancematrix', action ='store_true', help='(Optional) Calculate a distance matrix.')
    # parser_analysis.add_argument('--contactmap', action ='store_true', help='(Optional) Calculate a contact map.')
    # parser_analysis.add_argument('--ssassign', action ='store_true', help='(Optional) Calculate a secondary structure assignment matrix.')
    # parser_analysis.add_argument('--rg', action ='store_true', help='(Optional) Calculate a radius of gyration distribution.')
    # parser_analysis.add_argument('--dmax', action ='store_true', help='(Optional) Calculate a maximum distance distribution.')
    # parser_analysis.add_argument('--eed', action ='store_true', help='(Optional) Calculate an end-to-end distance distribution.')
    # parser_analysis.add_argument('--cmdist', type=dict[str,(str,str)], help='(Optional) Calculate a center of mass distance distribution for each pair of given MDAnalysis selection strings.', metavar='')
    # parser_analysis.add_argument('--colors', type=list[str], help='(Optional) List of color hexcodes to use, one for each analyzed trajectory.', metavar='')

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
    # parser_reweighting.add_argument('--rg', action ='store_true', help='(Optional) Calculate and compare uniform/reweighted radius of gyration distributions.')
    # parser_reweighting.add_argument('--dmax', action ='store_true', help='(Optional) Calculate and compare uniform/reweighted maximum distance distributions.')
    # parser_reweighting.add_argument('--eed', action ='store_true', help='(Optional) Calculate and compare uniform/reweighted end-to-end distance distributions.')
    # parser_reweighting.add_argument('--cmdist', type=dict[str,(str,str)], help='(Optional) Calculate and compare uniform/reweighted center of mass distance distributions for each pair of given MDAnalysis selection strings.', metavar='')

    # Now parse the remaining arguments with the full parser
    full_args = parser.parse_args([args.module] + remaining_args)

    # Handle the different modules based on the parsed arguments
    help_msg = '''
usage: ensemblify {generation,conversion,analysis,reweighting} [module options]

Command-line tool to access the modules of the Ensemblify Python library.

positional arguments:

    generation (g, gen)        Access the generation module.
    conversion (c, con)        Access the conversion module.
    analysis (a, ana)          Access the analysis module.
    reweighting (r, rew)       Access the reweighting module.
    clash_checking (cch)       Access the clash checking module.

'''
    if full_args.module in ['help','h']:
        print(help_msg)

    elif full_args.module in ['clash_checking', 'cch']:
        # from ensemblify.utils.clash_checking import check_steric_clashes
        # check_steric_clashes(ensemble_dir=full_args.ensembledir,
        #                      sampling_targets=full_args.samplingtargets,
        #                      input_structure=full_args.inputstructure)
        print(full_args.ensembledir)
        print(full_args.samplingtargets)
        print(full_args.inputstructure)
    elif full_args.module in ['generation', 'g', 'gen']:
        # from ensemblify import generate_ensemble
        # generate_ensemble(parameters_path=full_args.parameters)

        print(full_args.parameters)

    elif full_args.module in ['conversion', 'c', 'con']:
        # from ensemblify import ensemble2traj
        # ensemble2traj(job_name=full_args.jobname,
        #               ensemble_dir=full_args.ensembledir,
        #               trajectory_dir=full_args.trajectorydir,
        #               size=full_args.size)

        print(full_args.jobname)
        print(full_args.ensembledir)
        print(full_args.trajectorydir)
        print(full_args.size)

    elif full_args.module in ['analysis', 'a', 'ana']:
        # from ensemblify import analyze_trajectory
        # analyze_trajectory(trajectories=full_args.trajectory,
        #                    topologies=full_args.topology,
        #                    trajectory_ids=full_args.trajectoryid,
        #                    output_directory=full_args.outputdir,
        #                    ramachandran_data=full_args.ramadata,
        #                    distancematrices=full_args.distancematrix,
        #                    contactmatrices=full_args.contactmap,
        #                    ssassignments=full_args.ssassign,
        #                    rg=full_args.rg,
        #                    dmax=full_args.dmax,
        #                    eed=full_args.eed,
        #                    cm_dist=full_args.cmdist,
        #                    color_palette=full_args.colors)

        print(full_args.trajectory)
        print(full_args.topology)
        print(full_args.trajectoryid)
        print(full_args.outputdir)
        print(full_args.ramadata)
        print(full_args.distancematrix)
        print(full_args.contactmap)
        print(full_args.ssassign)
        print(full_args.rg)
        print(full_args.dmax)
        print(full_args.eed)
        print(full_args.cmdist)
        print(full_args.colors)

    elif full_args.module in ['reweighting', 'r', 'rew']:
        # from ensemblify import reweight_ensemble
        # reweight_ensemble(trajectory=full_args.trajectory,
        #                   topology=full_args.topology,
        #                   trajectory_id=full_args.trajectoryid,
        #                   exp_saxs_data=full_args.expdata,
        #                   output_dir=full_args.outputdir,
        #                   thetas=full_args.thetas,
        #                   compare_rg=full_args.rg,
        #                   compare_dmax=full_args.dmax,
        #                   compare_eed=full_args.eed,
        #                   compare_cmdist=full_args.cmdist)

        print(full_args.trajectory)
        print(full_args.topology)
        print(full_args.trajectoryid)
        print(full_args.expdata)
        print(full_args.outputdir)
        print(full_args.thetas)
        print(full_args.rg)
        print(full_args.dmax)
        print(full_args.eed)
        print(full_args.cmdist)
