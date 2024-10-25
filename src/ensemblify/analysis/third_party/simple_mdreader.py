"""
The code below was adapted from https://github.com/mnmelo/MDreader.
See the THIRD_PARTY_NOTICE_mdreader.txt file for more details.
"""

# IMPORTS
## Standard Library Imports
import datetime
import math
import multiprocessing
import sys

## Third Party Imports
import MDAnalysis as mda
import numpy as np

# CONSTANTS
# Default is to handle own errors, with a neat exit.
# Change to allow exceptions to reach the calling code.
RAISE_EXCEPTIONS = False

# FUNCTIONS
def parallel_launcher(
    rdr: mda.Universe,
    w_id: int):
    """ Helper function for the parallel execution of registered functions."""
    rdr.p_id = w_id
    return rdr._reader()

def raise_error(exc: Exception, msg: str):
    if RAISE_EXCEPTIONS:
        raise exc(msg)
    else:
        sys.exit(f'{exc.__name__}: {msg}')

# CLASSES
class Pool():
    """MDAnalysis and multiprocessing's map don't play along because of pickling.
    This solution seems to work fine.
    """
    def __init__(self, processes):
        self.nprocs = processes

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
                procs.append(multiprocessing.Process(target=self.fcaller,
                                                     args=(f,arglist.pop(0),num)))
                num += 1
                freeprocs -= 1
                procs[-1].start()
            # Execution halts here waiting for output after filling the procs.
            i, r = self.outqueue.get()
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

    def fcaller(self, f, args, num):
        self.outqueue.put((num, f(*args)))


class ThenNow:
    """Helper class to report analysis progress."""
    def __init__(self, oldval=None, newval=None):
        self.set(oldval, newval)

    def set(self, oldval, newval):
        self.old = oldval
        self.new = newval

    def fill(self, val):
        # Fill variant for the initial case where we have to assign both at initialization.
        self.set(val, val)

    def update(self, val, fill=False):
        if fill:
            self.fill(val)
        else:
            self.old = self.new
            self.new = val


# SimpleMDreader Class
class SimpleMDreader(mda.Universe):
    """An object class inheriting from MDAnalysis.Universe.
    
    Attributes:
        verbose (bool):
            whether to output messages to stdout reporting on analysis progress.
        nworkers (int):
            number of processor cores to use.
        nframes (int):
            size of loaded trajectory.
        _startframe (int):
            frame index from which to start analysis from. Defaults to 0.
        _endframe (int):
            frame index where to stop analysis. Defaults to self.nframes -1.
        totalframes (int):
            total number of frames to be analyzed.
        outstats (int):
            controls how often to report performance statistics. Defaults to 1, which leads to
            reporting statistics on every frame.
        statavg (int):
            controls over how many frames to accumulate performance statistics.
        loop_dtimes (np.ndarray):
            difference in times between current and previous self.statavg sized block of frames.
        loop_time (ThenNow):
            instance of custom helper class for reporting time statistics.
        framestr (str):
            formattable string to be printed at each time statistics report.
        p_num (int):
            actual number of processor cores that are being used. Defaults to all cores in machine.
        p_id (int):
            id of the current worker for reading and output purposes (to avoid terminal clobbering
            only p_id 0 will output). If messing with worker numbers (why would you do that?)
            beware to always set a different p_id per worker when iterating in parallel, otherwise
            you'll end up with repeated trajectory chunks.
        p_scale_dt (bool):
            controls whether the reported time per frame will be scaled by the number of workers,
            in order to provide an effective, albeit estimated, per-frame time.
        p_parms_set (bool):
            informs on whether the parallelization parameters have been setup or not.
        i_params_set (bool):
            informs on whether the iteration parameters have been setup or not.
    """
    def __init__(self,
        trajectory: str,
        topology: str,
        nworkers: int | None = None,
        outstats: int = 1,
        statavg: int = 100):
        """Initializes the SimpleMDReader instance based on the given parameters.
        
        Args:
            trajectory:
                path to trajectory file (.xtc).
            topology:
                path to topology file (.pdb).
            nworkers:
                number of processor cores to use during parallel calculations. If None, all
                cores in the machine are used.
            outstats:
                controls how often to report performance statistics.
            statavg:
                controls over how many frames to accumulate performance statistics.
        """
        self.verbose = True
        self.nworkers = nworkers
        mda.Universe.__init__(self, topology, trajectory)
        self.nframes = len(self.trajectory)
        self._startframe = 0
        self._endframe = self.nframes-1
        self.totalframes = int(np.rint(math.ceil(float(self._endframe - self._startframe+1))))

        # Stuff pertaining to progress output/parallelization
        self.outstats = outstats
        self.statavg = statavg
        self.loop_dtimes = np.empty(self.statavg, dtype=datetime.timedelta)
        self.loop_time = ThenNow()
        self.framestr = '{1:3.0%}  '
        self.p_num = nworkers
        self.p_id = 0
        self.p_scale_dt = True # set whether the per frame time considers the number of workers
        self.p_parms_set = False
        self.i_parms_set = False

    def p_fn(self):
        """The overridable function for parallel processing."""
        pass

    def _output_stats(self):
        """Keeps and outputs performance stats."""
        self.loop_time.update(datetime.datetime.now())
        
        if self.iterframe:  # No point in calculating delta times on iterframe 0
            self.loop_dtime = self.loop_time.new - self.loop_time.old
            self.loop_dtimes[(self.iterframe-1) % self.statavg] = self.loop_dtime
            
            # Output stats every outstat step or at the last frame
            if (not self.iterframe % self.outstats) or self.iterframe == self.i_totalframes - 1:
                avgframes = min(self.iterframe, self.statavg)
                self.loop_sumtime = sum(self.loop_dtimes[:avgframes], datetime.timedelta())
                
                # Compute ETA in seconds
                etaseconds = self.loop_sumtime.total_seconds() * (self.i_totalframes - self.iterframe) / avgframes
                eta = datetime.timedelta(seconds=etaseconds)
                
                if etaseconds > 300:
                    etastr = (datetime.datetime.now() + eta).strftime('Will end %Y-%m-%d at %H:%M:%S.')
                else:
                    etastr = f'Will end in {round(etaseconds)}s.'
                
                loop_dtime_s = self.loop_dtime.total_seconds()
                
                if self.p_scale_dt:
                    loop_dtime_s /= self.p_num

                progstr = self.framestr.format(self.snapshot.frame - 1, (self.iterframe + 1) / self.i_totalframes)
                
                sys.stderr.write(f'\033[K{progstr} ({loop_dtime_s:.4f} s/frame) \t{etastr}\r')
                
                if self.iterframe == self.i_totalframes - 1:
                    # Last frame, clean up
                    sys.stderr.write('\n')
                sys.stderr.flush()

    def _set_iterparms(self):
        # Because of parallelization lots of stuff become limited to the iteration scope.
        # defined a group of i_ variables just for that.
        self.i_unemployed = False

        # As-even-as-possible distribution of frames per workers, allowing the first one
        # to work more to compensate the lack of overlap.
        frames_per_worker = np.ones(self.p_num,dtype=int)*(self.totalframes//self.p_num)
        frames_per_worker[:self.totalframes%self.p_num] += 1
        self.i_startframe = int(self._startframe + np.sum(frames_per_worker[:self.p_id]))
        self.i_endframe = int(self.i_startframe + (frames_per_worker[self.p_id]-1))

        # Let's check for zero work
        if not frames_per_worker[self.p_id]:
            self.i_unemployed = True

        self.i_totalframes = int(np.rint(math.ceil((self.i_endframe-self.i_startframe+1))))
        self.i_parms_set = True

    def iterate(self):
        """Yields snapshots from the trajectory according to specified start and end boundaries.

        Calculations on AtomSelections will automagically reflect the new snapshot, without
        needing to refer to it specifically.

        SimpleMDreader.p_num, SimpleMDreader.p_id and SimpleMDreader.p_scale_dt are important
        here, read the class docstring for more information.
        """
        self.iterframe = 0

        # Let's always flush, in case the user likes to print stuff themselves.
        sys.stdout.flush()
        sys.stderr.flush()

        # The LOOP!
        for self.snapshot in self.trajectory[self.i_startframe:self.i_endframe+1]:
            if self.verbose and self.p_id==0:
                self._output_stats()
            yield self.snapshot
            self.iterframe += 1
        self.i_parms_set = False
        self.p_parms_set = False

    def _reader(self):
        """Applies self.p_fn for every trajectory frame. Parallelizable!"""
        # We need a brand new file descriptor per worker, otherwise we have a nice chaos.
        # This must be the first thing after entering parallel land.

        # XTC/TRR reader has this method, but not all formats...
        rdr = self.trajectory
        if hasattr(rdr, '_reopen'):
            rdr._reopen()
        else:
            raise_error(AttributeError, ('Don\'t know how to get a new file descriptor for '
                                         f'the {rdr.format} trajectory format. You\'ll have '
                                         'to skip parallelization.'))

        reslist = []
        if not self.i_parms_set:
            self._set_iterparms()
        if self.i_unemployed: # This little piggy stays home
            self.i_parms_set = False
            self.p_parms_set = False
            return reslist

        for _ in self.iterate():
            result = self.p_fn(*self.p_args, **self.p_kwargs)
            reslist.append(result)
        return reslist

    def do_in_parallel(self,
        fn,
        *args,
        **kwargs,
        ) -> list[float] | None:
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
        if self.nworkers is not None:
            self.p_num = self.nworkers
        else:
            self.p_num = multiprocessing.cpu_count()
        self.p_parms_set = True

        # Get args and kwargs for function to call
        self.p_args = args
        self.p_kwargs = kwargs

        pool = Pool(processes=self.p_num)
        res = pool.map(parallel_launcher, [(self, i) for i in range(self.p_num)])

        if self.p_id == 0:
            return [val for subl in res for val in subl]
