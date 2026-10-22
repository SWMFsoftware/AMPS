"""Wall-time limit of one fit: stop cleanly before a deadline (C15/P0; user decision 2026-10-22).

The fitter (``fit_cops``) stays unchanged. ``tools/gpu_fit.py`` replaces ``fit_cops.least_squares``,
the name ``fit_one`` calls, by ``WallClock.least_squares_with(original)``, which times every residual
and Jacobian evaluation. Before an evaluation it stops when

    now + (estimate of that evaluation) + reserve > deadline,
    reserve = save_seconds + estimate(residual) + estimate(Jacobian),

because after the optimizer ``fit_one`` evaluates the data residuals and the data Jacobian once more,
then the writers, the restart file and the figures run (``save_seconds``). The estimate of a kind is the
largest of its last ``HISTORY`` durations. The first evaluation of a kind has no estimate and runs, unless
not even ``save_seconds`` remain before the first residual (``WallTimeExhausted``: nothing to save).

On a stop, the optimizer result is the best evaluated state (lowest cost among all residual evaluations),
with ``status = STOP_STATUS`` (-2), ``success = False`` and the number of residual evaluations used.
``fit_one`` then reports it as usual; ``gpu_fit`` writes the restart file with ``segment.status =
"interrupted"`` and exits with ``EXIT_WALL_TIME`` (75). A later run continues from that restart file with
``--restart-from`` and the remaining evaluations. That is a new ``least_squares`` run (the trust region and
the Jacobian scaling start again), so **a segmented fit is not identical to one continuous run**.

Deadlines: ``--wall-time HH:MM:SS`` (from the start of the process), ``--wall-time slurm`` (the job's end:
``$SLURM_JOB_END_TIME``, else ``squeue -h -j $SLURM_JOB_ID -o %e``) or ``--deadline EPOCH`` (absolute, as
the pipeline runner passes it); a safety margin is subtracted. Default: no limit, and then the wrapper is
not installed at all.
"""

import datetime as dt
import os
import subprocess
import time

import numpy as np

EXIT_WALL_TIME = 75                 # EX_TEMPFAIL: stopped before the deadline, resumable
STOP_STATUS = -2                    # optimizer status of a fit stopped by the deadline
DEFAULT_MARGIN_SECONDS = 300.0      # safety margin before the deadline
DEFAULT_SAVE_SECONDS = 600.0        # reserved for the writers, restart file and figures (20M: ~6 min)
HISTORY = 5                         # evaluations per kind that define the estimate


class WallTimeExhausted(BaseException):
    """Not even the time to save one evaluation remains: stop without a result."""


class _Stop(Exception):
    """Internal: raised inside ``least_squares`` to end it at the best state so far."""


def parse_duration(text):
    """Seconds of ``[D-]HH:MM:SS``, ``MM:SS`` or plain seconds (Slurm's time formats)."""
    text = str(text).strip()
    days = 0
    if "-" in text:
        head, text = text.split("-", 1)
        days = int(head) if head.isdigit() else -1
    try:
        parts = [float(p) for p in text.split(":")]
    except ValueError:
        parts = []
    if not 1 <= len(parts) <= 3 or any(p < 0 for p in parts) or days < 0:
        raise ValueError("wall time must be [D-]HH:MM:SS, MM:SS or seconds: %r" % text)
    seconds = 0.0
    for p in parts:
        seconds = seconds * 60.0 + p
    return days * 86400.0 + seconds


def slurm_end_time(environ=None, run=subprocess.run):
    """End of the current Slurm job (epoch seconds)."""
    environ = os.environ if environ is None else environ
    value = environ.get("SLURM_JOB_END_TIME", "")
    if value.isdigit() and int(value) > 0:
        return float(value)
    job = environ.get("SLURM_JOB_ID")
    if not job:
        raise ValueError("--wall-time slurm outside a Slurm job (no SLURM_JOB_ID)")
    try:
        out = run(["squeue", "-h", "-j", job, "-o", "%e"], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                  universal_newlines=True, timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError) as exc:
        raise ValueError("cannot read the end time of job %s with squeue: %s" % (job, exc))
    try:
        end = dt.datetime.strptime(out.splitlines()[0].strip(), "%Y-%m-%dT%H:%M:%S")
    except (ValueError, IndexError):
        raise ValueError("job %s has no end time (squeue %%e gave %r)" % (job, out))
    return time.mktime(end.timetuple())            # squeue prints local time


def resolve_deadline(wall_time=None, deadline=None, margin=DEFAULT_MARGIN_SECONDS, start=None, environ=None):
    """Absolute stop time (epoch seconds, margin subtracted) or None (no limit)."""
    if deadline is not None:
        end = float(deadline)
    elif wall_time in (None, "", "none"):
        return None
    elif wall_time == "slurm":
        end = slurm_end_time(environ)
    else:
        end = (time.time() if start is None else start) + parse_duration(wall_time)
    return end - float(margin)


def utc(epoch):
    return None if epoch is None else dt.datetime.fromtimestamp(epoch, dt.timezone.utc).isoformat()


def scaled_optimality(x, g, lower, upper):
    """``least_squares``' first-order optimality for TRF with bounds: ||g * v||_inf (Coleman-Li v)."""
    v = np.ones_like(x)
    upper_side = (g < 0) & np.isfinite(upper)
    lower_side = (g > 0) & np.isfinite(lower)
    v[upper_side] = upper[upper_side] - x[upper_side]
    v[lower_side] = x[lower_side] - lower[lower_side]
    return float(np.max(np.abs(g * v))) if len(x) else 0.0


class WallClock(object):
    """Times the optimizer's evaluations and stops before ``deadline`` (see the module docstring)."""

    def __init__(self, deadline, save_seconds=DEFAULT_SAVE_SECONDS, clock=time.time):
        self.deadline = deadline
        self.save_seconds = float(save_seconds)
        self.clock = clock
        self.durations = {"fun": [], "jac": []}
        self.best = None                     # (cost, x, residual vector)
        self.last_jacobian = None            # (x, J) of the latest Jacobian evaluation
        self.nfev = self.njev = 0
        self.stopped = None

    def estimate(self, kind):
        recent = self.durations[kind][-HISTORY:]
        return max(recent) if recent else None

    def reserve(self):
        return self.save_seconds + (self.estimate("fun") or 0.0) + (self.estimate("jac") or 0.0)

    def check(self, kind):
        if self.deadline is None:
            return
        now = self.clock()
        estimate = self.estimate(kind)
        if estimate is None:
            if self.best is None and now + self.save_seconds > self.deadline:
                raise WallTimeExhausted("%.0f s before the deadline, %.0f s are needed to save a result"
                                        % (self.deadline - now, self.save_seconds))
            return                            # the first evaluation of a kind: just start it
        if now + estimate + self.reserve() > self.deadline:
            self.stopped = {"before": kind, "utc": utc(now), "seconds_left": self.deadline - now,
                            "estimate_seconds": estimate, "reserve_seconds": self.reserve()}
            raise _Stop()

    def timed(self, kind, function):
        def evaluate(x, *args, **kwargs):
            self.check(kind)
            t0 = self.clock()
            value = function(x, *args, **kwargs)
            self.durations[kind].append(self.clock() - t0)
            if kind == "fun":
                self.nfev += 1
                f = np.array(value, dtype=float, copy=True)
                cost = 0.5 * float(np.dot(f, f))
                if self.best is None or cost < self.best[0]:
                    self.best = (cost, np.array(x, dtype=float, copy=True), f)
            else:
                self.njev += 1
                self.last_jacobian = (np.array(x, dtype=float, copy=True), value)
            return value
        return evaluate

    def least_squares_with(self, original):
        """A drop-in ``least_squares`` that stops at the best state before the deadline."""
        def least_squares(fun, x0, jac="2-point", **kwargs):
            wrapped_jac = self.timed("jac", jac) if callable(jac) else jac
            try:
                return original(self.timed("fun", fun), x0, jac=wrapped_jac, **kwargs)
            except _Stop:
                return self.result(kwargs.get("bounds"))
        return least_squares

    def result(self, bounds=None):
        from scipy.optimize import OptimizeResult
        cost, x, f = self.best
        optimality = -1.0                     # not evaluated: no Jacobian at the best state
        if self.last_jacobian is not None and np.array_equal(self.last_jacobian[0], x):
            J = self.last_jacobian[1]
            g = np.asarray(J.T.dot(f), dtype=float).ravel()
            lower, upper = (np.broadcast_to(np.asarray(b, dtype=float), x.shape) for b in
                            (bounds if bounds is not None else (-np.inf, np.inf)))
            optimality = scaled_optimality(x, g, lower, upper)
        message = ("stopped before the wall-time deadline after %d residual evaluations; best evaluated "
                   "state (cost %.6e); resume with --restart-from" % (self.nfev, cost))
        return OptimizeResult(x=x, cost=cost, fun=f, jac=None, grad=None, optimality=optimality,
                              active_mask=np.zeros(len(x), dtype=int), nfev=self.nfev, njev=self.njev,
                              status=STOP_STATUS, message=message, success=False)

    def record(self):
        def stats(values):
            return ({"count": len(values), "mean": float(np.mean(values)), "max": float(np.max(values))}
                    if values else {"count": 0})
        return {"deadline_utc": utc(self.deadline), "save_seconds": self.save_seconds,
                "stopped": self.stopped, "nfev": self.nfev, "njev": self.njev,
                "best_cost": None if self.best is None else self.best[0],
                "evaluation_seconds": {k: stats(v) for k, v in self.durations.items()}}
