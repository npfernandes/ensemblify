"""python -m trace --trackcalls mdreader_callstack.py"""
# IMPORTS
# mdreader
from __future__ import division
import six
import sys
import argparse
import os
import numpy as np
import MDAnalysis
import math
import datetime
import multiprocessing

## Standard Library Imports
import os
import sys

## Third Party Imports
import MDAnalysis as mda
import numpy as np
import pandas as pd
import scipy

## Local Imports
#from ensemblify.analysis.third_party.mdreader_CHANGED import MDreader

class Tracker:
    level = 0

    def __init__(self, indent=2):
        self.indent = indent

    def __call__(self, fn):
        def wrapper(*args, **kwargs):
            print(' '*(self.indent * self.level) + '-' + fn.__name__,file=open('call_stack_output.txt','a',encoding='utf-8'))
            self.level += 1
            out = fn(*args, **kwargs)
            self.level -= 1
            return out
        return wrapper


track = Tracker()

# Globals ##############################################################
########################################################################
# Default is to handle own errors, with a neat exit. Change to allow
#  exceptions to reach the calling code.
raise_exceptions = False
_default_opts = {'s'    : 'topol.tpr',
                 'f'    : 'traj.xtc',
                 'o'    : 'data.xvg',
                 'b'    : 'from the beginning of the trajectory',
                 'e'    : 'until the end of the trajectory',
                 'skip' : 1,
                 'np'   : 0,
                 'v'    : 1}

# Helper functions and decorators ######################################
########################################################################

@track
def _with_defaults(defargs, clobber=True):
    """Decorator to set functions' default arguments

    The decorator takes a dictionary as an argument, which will
    supply default values for all the function arguments that match
    one of its keys.
    
    The decorator doesn't reorder function arguments. This means that, 
    as with def statements, after one argument has been assigned a default
    value all the following arguments must also be given default values.

    The decorator's clobber keyword defines whether existing default values
    are kept or clobbered by values in the passed dictionary.
    """
    def _fnmod(f):
        f_vars = f.__code__.co_varnames[:f.__code__.co_argcount]
        if f.__defaults__:
            ndefs = len(f.__defaults__)
            f_defaults = dict(zip(f_vars[-ndefs:], f.__defaults__)) 
        else:
            f_defaults = {}
        if clobber:
            f_defaults.update(defargs) 
        else:
            f_defaults, f_d = defargs.copy(), f_defaults
            f_defaults.update(f_d)
        new_defaults = []
        for var in f_vars:
            try:
                new_defaults.append(f_defaults[var])
            except KeyError:
                if new_defaults:
                    prev_arg = f_vars[f_vars.index(var)-1]
                    raise TypeError("While attempting to set defaults for the arguments of function "
                            "'{fname}' argument '{arg}' comes after optional argument '{prev_arg}' but was assigned "
                            "no default value. Either set a default value for '{arg}' or modify the base function "
                            "so that '{arg}' comes before any optional arguments.".format(fname=f.func_name, arg=var, prev_arg=prev_arg))
        f.__defaults__ = tuple(new_defaults)
        return f
    return _fnmod

@track
def _do_be_flags(val, default, asframenum):
    if val == default:
        return None
    else:
        if asframenum:
            val = int(val)
        else:
            val = float(val)
        #check_positive(val)
        return val

@track
def _parallel_launcher(rdr, w_id):
    """ Helper function for the parallel execution of registered functions.

    """
    rdr.p_id = w_id
    return rdr._reader()

@track
def raise_error(exc, msg):
    if raise_exceptions:
        raise exc(msg)
    else:
        sys.exit("{}: {}".format(exc.__name__, msg))

@track
def check_file(fname):
    if not os.path.exists(fname):
        raise_error(IOError, 'Can\'t find file %s' % (fname))
    if not os.access(fname, os.R_OK):
        raise_error(IOError, 'Permission denied to read file %s' % (fname))
    return fname

@track
def check_outfile(fname):
    dirname = os.path.dirname(fname)
    if not dirname:
        dirname = '.'
    if not os.access(dirname, os.W_OK):
        raise_error(IOError, 'Permission denied to write file %s' % (fname))
    return fname

@track
def check_positive(val, strict=False):
    if strict and val <= 0:
        raise_error(ValueError, "Argument '%r' must be > 0" % (val))
    elif val < 0:
        raise_error(ValueError, "Argument '%r' must be >= 0" % (val))

# Workaround for the lack of datetime.timedelta.total_seconds() in python<2.7
if hasattr(datetime.timedelta, "total_seconds"):
    dtime_seconds = datetime.timedelta.total_seconds
else:
    def dtime_seconds(dtime):
        return dtime.days*86400 + dtime.seconds + dtime.microseconds*1e-6

# Helper Classes #######################################################
########################################################################

class Pool():
    """ MDAnalysis and multiprocessing's map don't play along because of pickling. This solution seems to work fine.

    """

    @track
    def __init__(self, processes):
        self.nprocs = processes

    @track
    def map(self, f, argtuple):
        procs = []
        nargs = len(argtuple)
        result = [None]*nargs
        arglist = list(argtuple)
        self.outqueue = multiprocessing.Queue()
        freeprocs = self.nprocs
        num = 0
        got = 0
        while arglist:
            while arglist and freeprocs:
                procs.append(multiprocessing.Process(target=self.fcaller, args=((f, arglist.pop(0), num) )))
                num += 1
                freeprocs -= 1
                # procs[-1].daemon = True
                procs[-1].start()
            i, r = self.outqueue.get() # Execution halts here waiting for output after filling the procs.
            result[i] = r
            got += 1
            freeprocs += 1
        # Must wait for remaining procs, otherwise we'll miss their output.
        while got < nargs:
            i, r = self.outqueue.get()
            result[i] = r
            got += 1
        for proc in procs:
            proc.terminate()
        return result

    @track
    def fcaller(self, f, args, num):
        self.outqueue.put((num, f(*args)))


class ThenNow:

    @track
    def __init__(self, oldval=None, newval=None):
        self.set(oldval, newval)

    @track
    def set(self, oldval, newval):
        self.old = oldval
        self.new = newval

    @track
    def fill(self, val):
        # Fill variant for the initial case where we have to assign both at initialization.
        self.set(val, val)

    @track
    def update(self, val, fill=False):
        if fill:
            self.fill(val)
        else:
            self.old = self.new
            self.new = val


# MDreader Class #######################################################
########################################################################

# Effectivelly, MDreader will also inherit from either argparse.ArgumentParser
# or from DummyParser.
class MDreader(MDAnalysis.Universe):
    """An object class inheriting from both argparse.ArgumentParser and MDAnalysis.Universe. Should be initialized as for argparse.ArgumentParser, with additional named arguments:
    Argument 'arguments' should be passed the list of command line arguments; it defaults to sys.argv[1:], which is very likely what you'll want.
    Argument 'outstats' defines how often (framewise) to output frame statistics. Defaults to 1.
    Argument 'statavg' defines over how many frames to average statistics. Defaults to 100.
    Argument 'internal_argparse' lets the user choose whether they want to let MDreader take care of option handling. Defaults to True. If set to False, a set of default filenames and most other options (starttime, endtime, etc.) will be used. Check functions setargs and add_ndx on how to change these defaults, or directly modify the mdreader.opts object attributes.

    The command-line argument list will default to:
    
    usage: %prog% [-h] [-f [TRAJ [TRAJ ...]]] [-s TOPOL]
                  [-o OUT] [-b TIME/FRAME]
                  [-e TIME/FRAME] [-fmn] [-skip FRAMES]
                  [-np NPROCS] [-v LEVEL] [-n [INDEX]]

    optional arguments:
      -h, --help                    show this help message and exit
      -f [TRAJ [TRAJ ...]]  file    The trajectory to analyze. If multiple files
                                    they'll be analyzed concatenated. (default: traj.xtc)
      -s TOPOL              file    .tpr, .gro, or .pdb file with the same atom
                                    numbering as the trajectory. (default: topol.tpr)
      -o OUT                file    The main data output file. (default: data.xvg)
      -b TIME/FRAME         real    Time to begin analysis from. If -fmn is set,
                                    -b takes instead an int, as the starting frame number.
                                    (default: from the beginning of the trajectory)
      -e TIME/FRAME         real    Time to end analysis at. If -fmn is set, -e
                                    takes instead an int, as the end frame number.
                                    (default: until the end of the trajectory)
      -fmn                  bool    Whether to interpret -b and -e as frame
                                    numbers (0-based). (default: False)
      -skip FRAMES          int     Interval between frames when analyzing.
                                    (default: 1)
      -np NPROCS            int     Number of processes to parallelize over when
                                    iterating. 1 means serial iteration, and 0 uses the
                                    OS-reported number of cores. Ignored when using MPI,
                                    or when the script specifically sets the number of
                                    parallelization workers. (default: 0)
      -v LEVEL              enum    Verbosity level. 0:quiet, 1:progress 2:debug
                                    (default: 1)
      -n [INDEX]            file    Index file. Defaults to 'index.ndx' if the
                                    filename is not specified. If this flag is omitted
                                    altogether index information will be built from
                                    residue names. (default: None)

    where %prog% is the 'prog' argument as supplied during initialization, or sys.argv[0] if none is provided.
    
    After MDreader instantiation the values of the defaults to the arguments can be changed using MDreader.setargs() (also for setting/unsetting automatic file IO checking; see function documentation). If a 'ver' argument is passed to setargs it will be displayed as the program version, and a '-V'/'--version' option for that purpose will be automatically created.
    The arguments for an MDreader instance can also be added or overridden using the add_argument() method (see the argparse documentation).
    The iterate() method will iterate over the trajectory according to the supplied options, yielding frames as it goes. You'll probably want to use it as part of a for-loop header.
    argparse deprecates using the 'version' argument to __init__. If you need to set it, use the setargs method.
    
    """

    internal_argparse = True

    @track
    def __new__(cls, *args, **kwargs):
        bases = (cls,) + cls.__bases__
        newcls = type(cls.__name__, bases + (argparse.ArgumentParser,), {})
        return super(MDreader, newcls).__new__(newcls)

    @track
    def __init__(self, arguments=sys.argv[1:], outstats=1, statavg=100, *args, **kwargs):
        """ Sets up the MDreader object, but doesn't initialize most heavy stuff.

        Option parsing and topology/trajectory loading is left to be done on a need basis.
        keyword 'arguments' allows one to specify a custom list of CLI-like arguments.
        keyword 'outstats' controls how often to report performance statistics.
        keyword 'statavg' controls over how many frames to accumulate performance statistics.
        Finally, keyword 'internal_argparse' allows one to specify whether to use argparse for
        option parsing (set to True) or to use a DummyParser instead (set to False). In the
        latter case one must later supply all needed options by hand, via the setargs method.
        """
        self.arguments = arguments

        argparse.ArgumentParser.__init__(self, *args, **kwargs)
        self.check_files = True # Whether to check for readability and writability of input and output files.

        self.version = None
        self.setargs()
        self._parsed = False
        self.hasindex = False
        self._nframes = None
        # Stuff pertaining to progress output/parallelization
        self.parallel = False  # Whether to parallelize
        self.p_smp = False  # SMP parallelization (within the same machine, or virtual machine)
        self.p_mpi = False  # MPI parallelization
        self.outstats = outstats
        self.statavg = statavg
        self.loop_dtimes = np.empty(self.statavg, dtype=datetime.timedelta)
        self.loop_time = ThenNow()
        self.progress = None
        self.framestr = "{1:3.0%}  "
        self.p_mode = 'block'
        self.p_overlap = 0
        self.p_num = None
        self.p_id = 0
        self.p_scale_dt = True
        self.p_mpi_keep_workers_alive = False
        self.p_parms_set = False
        self.i_parms_set = False
        self._cdx_meta = False # Whether to also return time/box arrays when extracting coordinates.

        # Check whether we're running under MPI. Not failsafe, but the user should know better than to fudge with these env vars.
        mpivarlst = ['PMI_RANK', 'OMPI_COMM_WORLD_RANK', 'OMPI_MCA_ns_nds_vpid',
                     'PMI_ID', 'SLURM_PROCID', 'LAMRANK', 'MPI_RANKID',
                     'MP_CHILD', 'MP_RANK', 'MPIRUN_RANK']
        self.mpi = bool(sum([var in os.environ.keys() for var in mpivarlst]))

    # The overridable function for parallel processing.
    @track
    def p_fn(self):
        pass

    @track
    def __len__(self):
        return self.totalframes

    @track
    def __getattr__(self, name):
    # This is a fancy way of only doing important stuff when needed. Also, the user no longer needs to call do_parse, even to access MDAnalysis sub objects.
        if name in ['startframe','endframe','totalframes']:
            try:
                return getattr(self, "_"+name)
            except AttributeError:
                self._set_frameparms()
                return getattr(self, "_"+name)
        if name == "_anchor_uuid":
            raise AttributeError
        if not self._parsed:
            self.do_parse()
            return getattr(self, name)
        else:
            raise AttributeError("%r object has no attribute %r" % (self.__class__, name))

    @property
    def nframes(self):
        if self._nframes is None:
            self.ensure_parsed()

            self._nframes = len(self.trajectory)
            if self._nframes is None or self._nframes < 1:
                raise_error(IOError, 'No frames to be read.')

        return self._nframes

    @_with_defaults(_default_opts)
    def setargs(self, s, f, o, b, e, skip, np, v, version=None, check_files=None):
        """ This function allows the modification of the default parameters of the default arguments without having
            to go through the hassle of overriding the args in question. The arguments to this function will override the defaults
            of the corresponding options. These defaults are taken even when internal_argparse has been set to False.
            In the particular case of the 'o' and 'np' arguments one can pass 'None' to hide the option. Check the
            set_parallel_parms method on how to set a specific parallelization number.
            check_files (also accessible via MDreader_obj.check_files) controls whether checks are performed on the readability and
            writabilty of the input/output files defined here (default behavior is to check).
        """
        # Slightly hackish way to avoid code duplication
        parser = self #if self.internal_argparse else self._dummyopts
        # Note: MUST always use dest as a kwarg, to satisfy the DummyParser. Anything without 'dest' will be ignored by it (only relevant when the user sets internal_argparse to False)
        parser.add_argument('-f', metavar='TRAJ', dest="infile", default=f, nargs="*",
                help = 'file\tThe trajectory to analyze. If multiple files they\'ll be analyzed concatenated.')
        parser.add_argument('-s', metavar='TOPOL', dest="topol", default=s,
                help = 'file\t.tpr, .gro, or .pdb file with the same atom numbering as the trajectory.')
        if o is None:
            parser.add_argument('-o', metavar='OUT', dest='outfile', default='data.xvg',
                    help = argparse.SUPPRESS)
        else:
            parser.add_argument('-o', metavar='OUT', dest='outfile', default=o,
                    help = 'file\tThe main data output file.')
        parser.add_argument('-b', metavar='TIME/FRAME', dest='starttime', default=b,
                help = 'real\tTime to begin analysis from. If -fmn is set, -b takes '
                    'instead an int, as the starting frame number.')
        parser.add_argument('-e', metavar='TIME/FRAME', dest='endtime', default=e,
                help = 'real\tTime to end analysis at. If -fmn is set, -e takes '
                    'instead an int, as the end frame number.')
        parser.add_argument('-fmn',  action='store_true', dest='asframenum',
                help = 'bool\tWhether to interpret -b and -e as frame numbers (0-based).')
        parser.add_argument('-skip', metavar='FRAMES', type=int, dest='skip', default=skip,
                help = 'int \tInterval between frames when analyzing.')
        if np is None:
            parser.add_argument('-np', metavar='NPROCS', type=int, dest='parallel', default=_default_opts['np'],
                    help = argparse.SUPPRESS)
        else:
            parser.add_argument('-np', metavar='NPROCS', type=int, dest='parallel', default=np,
                    help = 'int \tNumber of processes to parallelize over when iterating. 1 means serial '
                    'iteration, and 0 uses the OS-reported number of cores. Ignored when using MPI, or when '
                    'the script specifically sets the number of parallelization workers.')
        parser.add_argument('-v', metavar='LEVEL', type=int, choices=[0,1,2], dest='verbose', default=v,
                help = 'enum\tVerbosity level. 0:quiet, 1:progress 2:debug')
        if version is not None:
            parser.add_argument('-V', '--version', action='version', version='%%(prog)s %s'%version,
                help = 'Prints the script version and exits.')
        if check_files is not None:
            self.check_files = check_files


    @track
    def ensure_parsed(self):
        if not self._parsed:
            self.do_parse()

    @track
    def do_parse(self):
        """ Parses the command-line arguments according to argparse and does some basic sanity checking on them. It also prepares some argument-dependent loop variables.
        If it hasn't been called so far, do_parse() will be called by the iterate() method, or when trying to access attributes that require it.
        Usually, you'll only want to call this function manually if you want to make sure at which point the arguments are read/parsed.

        """
        self.opts = self.parse_args(self.arguments)

        # We find the version string in the parser _actions. Somewhat fragile.
        for action in self._actions:
            if isinstance(action, argparse._VersionAction) and self.version is None:
                self.version = action.version

        # if self.opts.verbose and self.p_id == 0:
        #     sys.stderr.write("Loading...\n")

        ## Post option handling. outfile and parallel might be unset.
        if isinstance(self.opts.infile, six.string_types):
            self.opts.infile = [self.opts.infile,]
        if self.check_files:
            map(check_file, [self.opts.topol] + self.opts.infile)
            check_outfile(self.opts.outfile)
        check_positive(self.opts.skip, strict=True)
        check_positive(self.opts.parallel)

        # -b/-e flag handling:
        self.opts.starttime = _do_be_flags(self.opts.starttime, _default_opts['b'], self.opts.asframenum)
        self.opts.endtime = _do_be_flags(self.opts.endtime, _default_opts['e'], self.opts.asframenum)

        #if self.opts.endtime is not None and self.opts.endtime < self.opts.starttime:
        #    raise_error(ValueError, 'Specified end time/frame lower than start time/frame.')

        if not self.p_parms_set:
            self.set_parallel_parms(self.opts.parallel)
        MDAnalysis.Universe.__init__(self, self.opts.topol, *self.opts.infile)

        self.hastime = True
        if not hasattr(self.trajectory.ts, 'time') or self.trajectory.dt == 0.:
            # if not self.opts.asframenum and not self.p_id:
            #     sys.stderr.write("Trajectory has no time information. Will interpret limits as frame numbers.\n")
            self.hastime = False
            self.opts.asframenum = True

        self._parsed = True
        self._set_frameparms()

        if not self.p_id:  # we're root
            if self.hasindex:
                self._parse_ndx()


    @track
    def iterate(self, p=None):
        """Yields snapshots from the trajectory according to the specified start and end boundaries and skip.
        Calculations on AtomSelections will automagically reflect the new snapshot, without needing to refer to it specifically.
        Argument p sets the number of workers, overriding any already set. Note that MDreader is set to use all cores by default, so if you want serial iteration you must pass p=1.
        Other output and parallelization behavior will depend on a number of MDreader properties that are automatically set, but can be changed before invocation of iterate():
          MDreader.progress (default: None) can be one of 'frame', 'pct', 'both', 'empty', or None. It sets the output to frame numbers, %% progress, both, or nothing. If set to None behavior defaults to 'frame', or 'pct' when iterating in parallel block mode.
          MDreader.p_mode (default: 'block') sets either 'interleaved' or 'block' parallel iteration.
          When MDreader.p_mode=='block' MDreader.p_overlap (default: 0) sets how many frames blocks overlap, to allow multi frame analyses (say, an average) to pick up earlier on each block.
          MDreader.p_num (default: None) controls in how many blocks/segments to divide the iteration (the number of workers; will use all the processing cores if set to None) and MDreader.p_id is the id of the current worker for reading and output purposes (to avoid terminal clobbering only p_id 0 will output). 
            **
            If messing with worker numbers (why would you do that?) beware to always set a different p_id per worker when iterating in parallel, otherwise you'll end up with repeated trajectory chunks.
            **
          MDreader.p_scale_dt (default: True) controls whether the reported time per frame will be scaled by the number of workers, in order to provide an effective, albeit estimated, per-frame time.

        """
        verb = self.opts.verbose and (not self.parallel or self.p_id==0)
        # We're only outputting after each worker has picked up on the pre-averaging frames
        self.i_overlap = True
        self.iterframe = 0

        # Let's always flush, in case the user likes to print stuff themselves.
        sys.stdout.flush()
        sys.stderr.flush()

        # The LOOP!
        for self.snapshot in self.trajectory[self.i_startframe:self.i_endframe+1:self.i_skip]:
            if self.i_overlap and self.iterframe >= self.p_overlap:
                self.i_overlap = False # Done overlapping. Let the output begin!
            if verb:
                self._output_stats()
            yield self.snapshot
            self.iterframe += 1
        self.i_parms_set = False
        self.p_parms_set = False


    @track
    def _output_stats(self):
        """Keeps and outputs performance stats.
        """
        self.framestr = "{1:3.0%}  "
        self.loop_time.update(datetime.datetime.now())
        if self.iterframe: # No point in calculating delta times on iterframe 0
            self.loop_dtime = self.loop_time.new - self.loop_time.old
            self.loop_dtimes[(self.iterframe-1) % self.statavg] = self.loop_dtime
            # Output stats every outstat    s step or at the last frame.
            if (not self.iterframe % self.outstats) or self.iterframe == self.i_totalframes-1:
                avgframes = min(self.iterframe,self.statavg)
                self.loop_sumtime = self.loop_dtimes[:avgframes].sum()
                # No float*dt multiplication before python 3. Let's scale the comparing seconds and set the dt ourselves.
                etaseconds = dtime_seconds(self.loop_sumtime)*(self.i_totalframes-self.iterframe)/avgframes
                eta = datetime.timedelta(seconds=etaseconds)
                if etaseconds > 300:
                    etastr = (datetime.datetime.now()+eta).strftime("Will end %Y-%m-%d at %H:%M:%S.")
                else:
                    etastr = "Will end in %ds." % round(etaseconds)
                loop_dtime_s = dtime_seconds(self.loop_dtime)

                if self.p_scale_dt:
                    loop_dtime_s /= self.p_num

                if self.hastime:
                    progstr = self.framestr.format(self.snapshot.frame-1, (self.iterframe+1)/self.i_totalframes, self.snapshot.time)
                else:
                    progstr = self.framestr.format(self.snapshot.frame-1, (self.iterframe+1)/self.i_totalframes)

                sys.stderr.write("\033[K%s(%.4f s/frame) \t%s\r" % (progstr, loop_dtime_s, etastr))
                if self.iterframe == self.i_totalframes-1:
                    #Last frame. Clean up.
                    sys.stderr.write("\n")
                sys.stderr.flush()


    @track
    def do_in_parallel(self, fn, *args, **kwargs):
        """ Applies fn to every frame, taking care of parallelization details.
        
        Returns a list with the returned elements, in order.
        args and kwargs should be an iterable, resp. a dictionary, of arguments
            that will be passed (with the star, resp. double-star, operator) to
            fn. Default to the empty tuple and empty dict.
        parallel can be set to False to force serial behavior. Setting it to
            True forces default parallelization behavior, overriding previous
            settings of self.p_num.
        ret_type can be set to "last_per_worker" to specify that only the last
            frame result per worker be returned. This is useful when dealing
            with returned objects that are updated along the several frames.
        Refer to the documentation on MDreader.iterate() for information on
        which MDreader attributes to set to change default parallelization
        options.

        """
        # Set function to call
        self.p_fn = fn

        # Make sure we are using all the OS reported number of processors
        # Replace self.set_parallel_params(0)
        self.parallel = True
        self.p_mpi = False
        self.p_smp = True
        self.p_num = 1 # multiprocessing.cpu_count()
        self.p_parms_set = True

        # Get args and kwargs for function to call
        self.p_args = args
        self.p_kwargs = kwargs

        pool = Pool(processes=self.p_num)
        res = pool.map(_parallel_launcher, [(self, i) for i in range(self.p_num)])

        return [val for subl in res for val in subl]


    @track
    def _reader(self):
        """ Applies self.p_fn for every trajectory frame. Parallelizable!

        """
        # We need a brand new file descriptor per SMP worker, otherwise we have a nice chaos.
        # This must be the first thing after entering parallel land.
        self._reopen_traj()

        reslist = []
        if not self.i_parms_set:
            self._set_iterparms()
        if self.i_unemployed: # This little piggy stays home
            self.i_parms_set = False
            self.p_parms_set = False
            return reslist

        for frame in self.iterate():
            result = self.p_fn(*self.p_args, **self.p_kwargs)
            if not self.i_overlap:
                reslist.append(result)
        return reslist


    @track
    def _reopen_traj(self):
       # Let's make this generic and always loop over a list of formats. If it's the ChainReader then they all get in.
       rdrs = []
       if self.trajectory.format == "CHAIN":
           rdrs.extend(self.trajectory.readers)
       else:
           rdrs.append(self.trajectory)
       for rdr in rdrs:
           # XTC/TRR reader has this method, but not all formats...
           if hasattr(rdr, "_reopen"):
               rdr._reopen()
           elif hasattr(rdr, "dcdfile"):
               rdr.dcdfile.close()
               rdr.dcdfile = open(self.trajectory.filename, 'rb')
           else:
               raise_error(AttributeError, "Don't know how to get a new file descriptor for the %s trajectory format. You'll have to skip parallelization." % rdr.format)

    @track
    def _set_frameparms(self):
        if self.opts.starttime is None:
            self._startframe = 0
        elif self.opts.starttime < 0:
            self._startframe = self.nframes + self.opts.starttime
        else:
            self._startframe = self.opts.starttime
        #
        if self.opts.endtime is None:
            self._endframe = self.nframes-1
        elif self.opts.endtime < 0:
            self._endframe = self.nframes + self.opts.starttime
        else:
            self._endframe = min(self.nframes-1, self.opts.endtime)

        if self._startframe >= self.nframes:
            raise_error(ValueError, "You requested to start at frame %d but the trajectory only has %d frames." % (self.opts.starttime, self.nframes))
        if self._endframe < self._startframe:
            raise_error(ValueError, 'Specified end frame lower than start frame.')
        if self._endframe < 0 or self._startframe < 0:
            raise_error(ValueError, 'Resulting start/end frame lower than 0.')

        self._totalframes = int(np.rint(math.ceil(float(self._endframe - self._startframe+1)/self.opts.skip)))

    @track
    def _set_iterparms(self):
        # Because of parallelization lots of stuff become limited to the iteration scope.
        # defined a group of i_ variables just for that.
        self.i_unemployed = False
        #if self.p_num < 2 and self.p_smp:
        #    raise ValueError("Parallel iteration requested, but only one worker (MDreader.p_num) sent to work.")

        # Assuming self.p_mode = 'block'

        #self.i_startframe and self.i_endframe must be specifically cast as
        # ints because of an unpythonic type check in MDAnalysis that misses numpy ints.
        # As-even-as-possible distribution of frames per workers, allowing the first one to work more to compensate the lack of overlap.
        frames_per_worker = np.ones(self.p_num,dtype=int)*((self.totalframes-self.p_overlap)//self.p_num)
        frames_per_worker[:(self.totalframes-self.p_overlap)%self.p_num] += 1
        frames_per_worker[0] += self.p_overlap # Add extra overlap frames to the first worker.
        self.i_skip = self.opts.skip
        self.i_startframe = int(self.startframe + np.sum(frames_per_worker[:self.p_id])*self.i_skip)
        self.i_endframe = int(self.i_startframe + (frames_per_worker[self.p_id]-1)*self.i_skip)
        # And now we subtract the overlap from the startframe, except for worker 0
        if self.p_id:
            self.i_startframe -= self.p_overlap*self.i_skip

        # Let's check for zero work
        if not frames_per_worker[self.p_id]:
            self.i_unemployed = True

        self.i_totalframes = int(np.rint(math.ceil((self.i_endframe-self.i_startframe+1)/self.i_skip)))
        self.i_parms_set = True


# FUNCTIONS
@track
def calc_rg(u: mda.Universe) -> float:
    """Calculate the radius of gyration of the current frame.
    
    Args:
        u:
            Universe pointing to the current frame.
    
    Returns:
        rg:
            radius of gyration of the protein in the current frame.
    """
    protein = u.select_atoms('protein')
    rg = protein.radius_of_gyration()
    return rg

@track
def calc_eed(u: mda.Universe) -> float:
    """Calculate the distance from the N to the C terminal in the current frame.
    
    Args:
        u:
            Universe pointing to the current frame.
    
    Returns:
        eed:
            end-to-end distance of the protein in the current frame.

    """
    nterm = u.select_atoms('protein and name N')[0]
    cterm = u.select_atoms('protein and name C')[-1]
    eed = np.linalg.norm(cterm.position - nterm.position)
    return eed

@track
def calc_dmax(u: mda.Universe) -> float:
    """Calculate the maximum of the distances between any two alpha carbons in the current frame.

    Args:
        u:
            Universe pointing to the current frame.
    
    Returns:
        dmax:
            Maximum of the distances between any two alpha carbons of the protein in the current
            frame.

    """
    ca_selection = u.select_atoms('protein and name CA')
    ca_coordinates = ca_selection.positions #expose numpy array of coords
    distance_matrix_pool = scipy.spatial.distance.cdist(ca_coordinates, ca_coordinates)
    maximum_distance_pool = distance_matrix_pool.max()
    dmax = np.linalg.norm(maximum_distance_pool)
    return dmax

@track
def calc_cm_dist(
    u: mda.Universe,
    sel1: str,
    sel2: str,
    ) -> float:
    """Calculate the distance between the center of mass of two atom selections
    in the current frame.

    Args:
        u:
            Universe pointing to the current frame.
        sel1:
            MDAnalysis selection string for selecting an AtomGroup whose center of mass will be
            calculated.
        sel2:
            MDAnalysis selection string for selecting an AtomGroup whose center of mass will be
            calculated.
    
    Returns:
        cm_dist:
            Center of mass distance between AtomGroups selected by sel1 and sel2.

    """
    cm1 = u.select_atoms(sel1).center_of_mass()
    cm2 = u.select_atoms(sel2).center_of_mass()
    cm_dist = np.linalg.norm(cm1 - cm2)
    return cm_dist

@track
def calculate_metrics_data(
    trajectory: str,
    topology: str,
    rg: bool = True,
    dmax: bool = True,
    eed: bool = True,
    cm_dist: dict[str,tuple[str,str]] | None = None,
    output_path: str | None = None,
    ) -> pd.DataFrame:
    """Calculate structural metrics for each frame of a trajectory.

    Args:
        trajectory:
            path to .xtc trajectory file.
        topology:
            path to .pdb topology file.
        rg:
            whether to calculate the radius of gyration of the protein.
        dmax:
            whether to calculate the maximum distance between any two alpha carbons in the protein.
        eed:
            whether to calculate the distance from the N to C terminal of the protein.
        cm_dist:
            mapping of identifiers to tuples with two selection strings for creating MDAnalysis
            AtomGroups, whose center mass distance will be calculated. For example:
                {'inter_domain' : ('resid 1:30', 'resid 110:140')}
            If None, no center mass distances are calculated.
            See https://userguide.mdanalysis.org/stable/selections.html for more information about
            MDAnalysis selections.
        output_path:
            path to output .csv file or output directory. If directory, written file is named
            'structural_metrics.csv'.

    Returns:
        traj_analysis:
            DataFrame where columns are the desired structural metrics and rows are the frames
            of the trajectory.
    """
    # Initialize MDReader
    u = MDreader(f'-f {trajectory} -s {topology}'.split())

    # Calculate trajectory metrics
    results = []
    if rg:
        print('Calculating rg...')
        rgs = np.array(u.do_in_parallel(calc_rg,u))
        results.append(('rg',rgs))
    if eed:
        print('Calculating eed...')
        eeds = np.array(u.do_in_parallel(calc_eed,u))
        results.append(('eed',eeds))
    if dmax:
        print('Calculating dmax...')
        dmaxs = np.array(u.do_in_parallel(calc_dmax,u))
        results.append(('dmax',dmaxs))
    if cm_dist:
        for cm_dist_id in cm_dist:
            print(f'Calculating {cm_dist_id}...')
            cm_dists = np.array(u.do_in_parallel(calc_cm_dist,
                                                 u,
                                                 cm_dist[cm_dist_id][0],
                                                 cm_dist[cm_dist_id][1]))
            results.append((cm_dist_id,cm_dists))

    # Extract column names and values
    column_ids = []
    values = []
    for metric_id,metric_values in results:
        column_ids.append(metric_id)
        values.append(metric_values)

    # Create trajectory analysis DataFrame
    metrics_array = np.dstack(tuple(values))
    traj_analysis = pd.DataFrame(metrics_array.reshape(-1,
                                                       metrics_array.shape[-1]),
                                                       columns=column_ids)

    if output_path is not None:
        # Save structural metrics
        if os.path.isdir(output_path):
            traj_analysis.to_csv(os.path.join(output_path,'structural_metrics.csv'))
        elif output_path.endswith('.csv'):
            traj_analysis.to_csv(output_path)
        else:
            print(('Structural metrics DataFrame was not saved to disk, '
                   'output path must be a directory or .csv filepath!'))

    return traj_analysis

#sys.setprofile(tracefunc)
calculate_metrics_data(trajectory='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/Hst5/Hst5_trajectory.xtc',
                       topology='/home/tiagogomes/Desktop/projects/nuno_fernandes/Ensembles_Without_AlphaFold/TRAJECTORIES/Hst5/Hst5_top.pdb',
                       output_path=os.getcwd(),
                       #dmax=False,
                       #eed=False
                       )