# Submitting parallel EGSnrc jobs

EGSnrc simulations are "embarrassingly parallel": every history is
statistically independent, so a simulation can be split into `N` separate
jobs, each running `1/N` of the histories with an independent random number
sequence, and the partial results combined at the end. All of the scripts
that drive this live in `$HEN_HOUSE/scripts`.

There are two systems in the distribution:

- **`egs-parallel`** — the current, actively maintained front-end
  (2020 onward). Recommended for new work.
- **`exb` / `run_user_code_batch`** — the original (2003 onward) front-end,
  still installed and aliased by default. Useful on sites that already have
  it configured, or for batch systems not covered by `egs-parallel`.

Both systems rely on the **built-in parallel-processing support of the
EGSnrc applications themselves** (selected with the `-P`, `-j`, `-f` and
`-s` command-line flags) — they are merely convenient wrappers that launch
multiple copies of the application with the right flags on a given batch
system. You can always invoke an application with these flags by hand
without using either wrapper.

## 1. How EGSnrc parallel runs work

An EGSnrc (egs++) application understands these command-line options:

| Option | Meaning |
|---|---|
| `-i input_file` | input file (without `.egsinp`) |
| `-p pegs_file` | PEGS4 data file (omit, or use `pegsless`, for pegsless mode) |
| `-b` | batch mode: redirect output to `output_file.egslog` |
| `-P N` | this run is job `i` of `N` parallel jobs |
| `-j i` | the index `i` (1..N) of this job |
| `-f first` | the index of the "first" job, the one responsible for creating the lock file (normally `1`) |
| `-s` / `--simple-run` | force a *simple* run-control object (RCO) instead of the default job-control-file (JCF) RCO |
| `-u` / `--urc` | force a *uniform* run-control object (RCO); takes precedence over `-s` if both are given |

When `-P`/`-j` are given, the application picks one of three run-control
objects (RCOs) to coordinate the parallel jobs:

| RCO | Selected by | Lock file? | History split | Combination |
|---|---|---|---|---|
| **balanced** (`EGS_JCFControl`, the default) | nothing special needed | yes, `basename.lock` | dynamic: jobs pull chunks from the lock file as they finish previous ones, so faster workers do more work | the last job to finish does it |
| **simple** (`EGS_RunControl`) | `-s` / `--simple-run`, or `rco type = simple` | no | static: each job gets a fixed, equal `1/N` share | none automatic — combine by hand (`calculation = combine`, see below) |
| **uniform** (`EGS_UniformRunControl`) | `-u` / `--urc`, or `rco type = uniform` | no | static: each job gets a fixed, equal `1/N` share, same as simple | automatic: one or more "watcher" jobs poll for the other jobs' `*_w<i>.egsdat` files and combine incrementally as they appear |

The **balanced/JCF** RCO is the historical default for parallel runs: it
gives the best load balancing (the input's `nchunk` parameter controls how
many chunks the histories are split into, default `nchunk*N`), and jobs do
not need to start or finish at the same time. Its downside is that it
depends on working file locking on the `EGS_HOME` filesystem — on some
shared/network filesystems (or when jobs are not guaranteed to start in
order, e.g. some cluster schedulers) the lock file can be unreliable.

The **simple** RCO sidesteps file locking entirely by giving every job a
fixed, equal share of the histories up front, but nothing combines the
results automatically — you have to run a separate `calculation = combine`
pass yourself once all jobs are done.

The **uniform** RCO (added to remove the dependency on `*.lock` files,
see [PR #588](https://github.com/nrc-cnrc/EGSnrc/pull/588)) also splits
histories into a fixed, equal `1/N` share per job like the simple RCO, but
restores automatic combination without needing a lock file: by default the
*last* job (`ipar == ifirst + npar - 1`) acts as a "watcher" — once it
finishes its own histories it polls, for up to `number of intervals` ×
`interval wait time` (default 5 × 1000 ms = 5 s total), for the other jobs'
`<final_output_file>_w<i>.egsdat` files, combining incrementally as they
appear. Because every job runs the same fixed number of histories, the
uniform RCO is only as fast as the slowest job (no dynamic load balancing),
which makes it best suited to homogeneous compute environments (identical
hardware/software per job) — hence `EGS_UniformRunControl`'s description in
the source as a RCO "for homogeneous computing environments (HCE)".

Uniform-RCO behaviour can be tuned via the `run control` input block:

| Key | Meaning | Default |
|---|---|---|
| `rco type = simple\|uniform\|balanced` | select the RCO when no `-s`/`-u` flag is given | `balanced` |
| `watcher jobs = j1,...,jn` | which job indices act as watchers (instead of just the last job) | last job only |
| `interval wait time = ms` | how long a watcher sleeps between completion checks | `1000` |
| `number of intervals = n` | how many times a watcher checks before giving up | `5` |
| `check jobs completed = yes\|no` | whether watchers actively poll for `*.egsdat` files at all | `yes` |

> **Note:** the in-source documentation comment for `check jobs completed`
> states its default is `no`, but reading the actual implementation
> (`egs_run_control.cpp`) shows it defaults to `yes` (checking enabled) when
> the key is omitted — the comment in the header appears to be stale. Set it
> explicitly if the behaviour matters to you.

If neither `-s`/`-u`/`--simple-run`/`--urc` is given on the command line and
the application's input has no `run control` block at all, parallel runs
(`-P` > 1) default to the balanced/JCF RCO and serial runs default to
simple — i.e. existing input files behave exactly as before these RCOs were
added.

To switch an `egs-parallel` or `exb` run from the default lock-file-based
behaviour to the uniform RCO, simply add `-u` (or `--urc`) to the
application invocation passed via `--command` (for `egs-parallel`) or the
trailing arguments (for `exb`/`egs-jobsub-xargs`/`egs-slurm`) — no other
changes are required, since `-P`/`-j`/`-f` are handled identically by all
three RCOs.

Two state files appear in the application directory while jobs are running:

- `basename.lock` — the job-control file (chunk bookkeeping). Not created
  when using the simple or uniform RCO.
- `basename.egsjob` — written by `egs-parallel`'s batch scripts; records
  `BEGIN`/`END` markers with host and PID, used to detect whether a
  simulation is still running or already finished.

These, plus `basename.egsparallel` (the log of the submission scripts
themselves), are removed/recreated automatically when you resubmit, and
cleaned up by `egs-parallel-clean` (see §4). Note that `egs-parallel`'s
batch scripts also stagger job startup and wait for `basename.egsjob` to
appear before launching jobs 2..N — this staggering exists to avoid a race
on `basename.lock` and is unnecessary (though harmless) for the uniform
RCO, since uniform jobs do not share any state until the watcher phase.

## 2. The modern way: `egs-parallel`

### 2.1 Overview

```
$HEN_HOUSE/scripts/bin/egs-parallel [options] --command 'command'
```

`egs-parallel` parses the command, figures out the application, input file
and "first job" index, then `exec`s a batch-specific helper script,
`$HEN_HOUSE/scripts/egs-parallel-<batch>`, which actually launches the `N`
copies of the application. `egs-parallel-clean` cleans up afterwards.

Make sure `$HEN_HOUSE/scripts/bin` is on your `PATH` (the
`egsnrc_bashrc_additions` / `egsnrc_cshrc_additions` / `egsnrc_fishrc_additions`
files sourced into your shell rc do this for you).

### 2.2 Options

```
egs-parallel [options] --command 'command'

  -h | --help         show help
  -b | --batch        batch system to use ("cpu" by default)
  -d | --delay        delay in seconds between individual jobs
  -q | --queue        scheduler queue ("long" by default)
  -n | --nthread      number of threads/jobs ("8" by default)
  -o | --option       option(s) to pass to the job scheduler, in quotes
  -f | --force        force run even if a lock or .egsjob file is present
  -v | --verbose      echo detailed egs-parallel log messages to terminal
  -c | --command      command to run, given in quotes
```

`--command` must contain a valid invocation of an EGSnrc executable with at
least `-i input_file`; `egs-parallel` adds `-P`, `-j` and `-f` itself, so do
not include those. The first job index defaults to `1`; pass `-f` inside
`--command` to override it (e.g. when resuming a partially-run set of
jobs).

Available batch systems (`-b`/`--batch`), each implemented as
`egs-parallel-<name>`:

| Batch | Script | Description |
|---|---|---|
| `cpu` | `egs-parallel-cpu` | Run `N` jobs as local background processes (one multicore machine, no scheduler). |
| `pbs` | `egs-parallel-pbs` | Submit `N` individual PBS jobs, one `qsub` per job. |
| `jobarr` | `egs-parallel-jobarr` | Submit one combined PBS **job array** (`#PBS -J`) covering all `N` jobs. |
| `pdsh` | `egs-parallel-pdsh` | Submit one PBS job that reserves `N` nodes and fans out tasks with `pdsh`/`ssh` (calls `egs-parallel-dshtask` on each node). |
| `pbsdsh` | `egs-parallel-pbsdsh` | Submit one PBS job that reserves `N` cores and fans out tasks with `pbsdsh` (calls `egs-parallel-dshtask` on each task). |
| `slurm` | `egs-parallel-slurm` | Submit `N` individual Slurm jobs, one `sbatch` per job. Local addition, modeled directly on `egs-parallel-pbs`. |

`egs-parallel-slurm` mirrors `egs-parallel-pbs` job-for-job: it submits one
`sbatch --parsable -p $queue $scheduler_options` per job (instead of one
`qsub` per job), each job body running the command and writing
`BEGIN`/`END` markers into `basename.egsjob` (job 1 only), with the same
staggered-start/lock-file logic as the PBS backend. Two notable differences
from `egs-parallel-pbs`:

- it uses `sbatch --parsable` to get a bare job ID back (no PBS-version
  branching is needed, since this output format is stable across Slurm
  versions);
- the job script does `cd $SLURM_SUBMIT_DIR` and exports the environment
  with `#SBATCH --export=ALL`, instead of PBS's `$PBS_O_WORKDIR` and
  `#PBS -v HEN_HOUSE,EGS_HOME,EGS_CONFIG`.

Pass extra `sbatch` flags (e.g. account, QOS, time limit, memory) with
`-o`/`--option`, the same way you would pass extra `qsub` flags to the PBS
backend:

```sh
egs-parallel --batch slurm -q standard -n 50 \
    -o '--account=myproject --qos=low --time=02:00:00' \
    --command 'egs_chamber -i slab -p 521icru'
```

### 2.3 Examples

Run 12 local jobs on a multicore workstation:

```sh
cd $EGS_HOME/egs_chamber
egs-parallel --batch cpu -n 12 -v --command 'egs_chamber -i slab -p 521icru'
```

Submit 50 individual PBS jobs to the `short` queue, with a 2 second delay
between submissions, passing extra `qsub` options:

```sh
egs-parallel --batch pbs -q short -n 50 -d 2 \
    -o '-l walltime=02:00:00' \
    --command 'egs_chamber -i slab -p 521icru'
```

Submit the same as a single PBS job array:

```sh
egs-parallel --batch jobarr -q short -n 50 \
    --command 'egs_chamber -i slab -p 521icru'
```

Reserve 20 nodes with `pdsh` under one PBS job:

```sh
egs-parallel --batch pdsh -q long -n 20 \
    --command 'dosxyznrc -i myphantom -p 521icru'
```

By default, `egs-parallel` refuses to run if `basename.lock` or
`basename.egsjob` already exist (it assumes a simulation is in progress or
finished). Add `-f`/`--force` to override, e.g. to resubmit after a crash.

### 2.4 What happens under the hood

1. `egs-parallel` validates the command, locates `basename.egsinp` under
   `$EGS_HOME/<app>`, and writes its own log to
   `$EGS_HOME/<app>/basename.egsparallel`.
2. It `exec`s `egs-parallel-<batch> queue nthread delay first basename
   'command' 'scheduler_options' verbosity`.
3. That script removes any stale `.egsjob`/`.lock` files, then launches `N`
   copies of `command -b -P N -j <1..N> -f <first>`:
   - `cpu`: as local background processes via `&`/`wait`.
   - `pbs`: one `qsub` call per job, each job body running the command and
     writing `BEGIN`/`END` markers into `basename.egsjob` (job 1 only).
   - `jobarr`: one `qsub` call with `#PBS -J first-nthread:1`; each array
     task derives its own job index from `$PBS_ARRAY_INDEX`.
   - `pdsh`/`pbsdsh`: one `qsub` call that reserves `N` nodes/cores, then
     fans out via `pdsh`/`pbsdsh` to `egs-parallel-dshtask`, which elects a
     "manager" task to hand out job indices to the others.
4. Every job, after the first, waits for `basename.egsjob` to appear
   (created by job 1) before starting, with small staggered delays (plus
   any user `--delay`) to avoid a thundering-herd on the lock file. Jobs
   also poll `basename.egsjob` for an `END` marker and quit immediately if
   the simulation has already finished (useful when resubmitting extra
   jobs into an already-completed run).
5. All activity is logged to `basename.egsparallel`. Run with `-v` to also
   see it on the terminal.

## 3. The legacy way: `exb` / `run_user_code_batch`

`exb` is a shell alias (defined in `egsnrc_bashrc_additions`,
`egsnrc_cshrc_additions`, `egsnrc_fishrc_additions`) for
`$HEN_HOUSE/scripts/run_user_code_batch`. It predates `egs-parallel` and is
still the way to reach batch systems that don't have an
`egs-parallel-<batch>` script (e.g. NQS; Slurm is now covered by
`egs-parallel --batch slurm`, see §2.2 above, but `exb batch=slurm` remains
available too).

### 3.1 Usage

```
exb user_code input_file pegs_file [noopt] [config=xxx]
    [short|long|user1|user2|user3] [batch=batch_system] [p=N]
    [start=S] [stop=E] [sleep=secs] [simple] [fresh=no]
```

- `user_code` — application name (e.g. `dosxyznrc`, `egs_chamber`).
- `input_file` — name without `.egsinp` (use `""` if none).
- `pegs_file` — PEGS4 data file name, or `pegsless`.
- `short|long|user1|user2|user3` — selects one of the queues defined in the
  batch-options file (see §3.2).
- `batch=xxx` — batch system to use (selects
  `$HEN_HOUSE/scripts/batch_options.xxx`); defaults to the environment
  variable `EGS_BATCH_SYSTEM`, or `at` if that isn't set either.
- `p=N` — submit `N` parallel jobs (adds `-P N -j <job> -f <first>` to the
  command for each job submitted).
- `start=S stop=E` — submit jobs numbered `S+1..E` instead of `1..N`
  (used to add extra jobs to, or resume part of, an existing run).
- `simple` — adds `-s` (force simple RCO).
- `fresh=no` — omit `-f` (don't treat this submission as starting a fresh
  job-control file; use when adding more jobs to an already-running set).
- `sleep=secs` — extra delay between individual `pbsdsh` task starts.

### 3.2 Batch system definitions: `batch_options.*`

`run_user_code_batch` sources `$HEN_HOUSE/scripts/batch_options.$batch`,
which defines:

- `batch_command` — the submission command (`qsub`, `sbatch`, `at`, ...).
- `generic_bo` — batch options independent of job name.
- `output_bo` / `rname_bo` — flags to redirect output and name the
  request, respectively.
- `batch_sleep_time` — delay between successive submissions (needed by
  some schedulers, e.g. PBS, to avoid getting confused by rapid submission).
- `batch_mnl` — maximum job-name length, truncated automatically.
- `short_queue` / `long_queue` / `user1_queue` / `user2_queue` /
  `user3_queue` — the actual queue selectors passed to `batch_command`.

Definitions shipped with the distribution:

| File | Scheduler |
|---|---|
| `batch_options.at` | Unix `at` (the default if `EGS_BATCH_SYSTEM` is unset) |
| `batch_options.nqs` | NQS |
| `batch_options.pbs` | PBS (`qsub`) |
| `batch_options.pbsdsh` | PBS distributed shell (`pbsdsh`) |
| `batch_options.slurm` | Slurm (`sbatch`) |
| `batch_options.keg` | site-specific (KEG) example |

To support a new scheduler, copy one of these to `batch_options.xxx`, adjust
the variables, and either `export EGS_BATCH_SYSTEM=xxx` or pass `batch=xxx`
to `exb`.

### 3.3 Examples

Submit a single run to the default (`at`) system:

```sh
exb dosxyznrc myphantom 521icru
```

Submit 20 parallel PBS jobs to the `short` queue:

```sh
exb egs_chamber slab 521icru short batch=pbs p=20
```

Submit to Slurm:

```sh
export EGS_BATCH_SYSTEM=slurm
exb egs_chamber slab 521icru long p=20
```

Reserve nodes via PBS's distributed shell (`pbsdsh`) — one PBS job that
internally runs `run_pbsdsh_task` on each reserved node:

```sh
exb egs_chamber slab 521icru batch=pbsdsh p=20
```

Add 5 more jobs (`21`..`25`) to an already-running set of 20, without
recreating the lock file:

```sh
exb egs_chamber slab 521icru batch=pbs p=5 start=20 stop=25 fresh=no
```

### 3.4 `egs-jobsub-xargs`

`$HEN_HOUSE/scripts/egs-jobsub-xargs` is a site-specific (NRC) variant for
Grid Engine, which batches several job indices into one scheduler
submission and fans them out locally with `xargs -P`:

```
egs-jobsub-xargs app inp pegs|pegsless [p=N] [config=xxx] [cpus=N] ...
```

It relies on a local `jobsub`/Grid-Engine setup and NRC-specific defaults
(project name, email, image, resource requests) hard-coded near the top of
the script — review/edit those before using it elsewhere.

### 3.5 `egs-slurm` (local, not part of the distribution)

`$HEN_HOUSE/scripts/egs-slurm` is a local NRC GPSC-cluster adaptation of
`egs-jobsub-xargs` for Slurm: instead of one job per `-j` index, it batches
`cpus` job indices (one node's worth) per `sbatch` submission and fans them
out inside that allocation with `xargs -P${cpus}`:

```
egs-slurm app inp pegs|pegsless [p=N] [config=xxx] [cpus=N] [cluster=xxx]
          [qos=xxx] [image=xxx] ...
```

Like `egs-jobsub-xargs`, it has NRC-specific defaults hard-coded near the
top (project, email, GPSC cluster/image/QOS, per-cluster core counts) —
review/edit those before using it elsewhere. Note a few rough edges found
while reviewing it for this document, worth fixing if you plan to keep
using it:

- the `cpus=N` command-line override is parsed into `cpus_u` but never
  actually used — `cpus` is unconditionally overwritten afterwards based on
  `$cluster`;
- the `rmegsdat` cleanup (`rm -f ...${inpf}_w*.egsdat`) runs before `inpf`
  is computed, so it always expands to an empty pattern (same pre-existing
  bug as in `egs-jobsub-xargs`);
- `executable="${egs_home}bin/$my_machine/$app"` is missing a `/` between
  `$egs_home` and `bin` (also inherited from `egs-jobsub-xargs`);
- the file header comment still says "Grid Engine job scheduler using
  xargs", left over from copying `egs-jobsub-xargs`.

For a non-NRC-specific, `egs-parallel`-integrated way to submit to Slurm,
see `egs-parallel-slurm` (§2.2) instead — it follows the same `-P`/`-j`/`-f`
one-job-per-submission model as the rest of `egs-parallel`'s backends,
rather than `egs-slurm`'s xargs-batching-per-node model.

### 3.6 `pprocess` (very old, mortran-era)

`$HEN_HOUSE/scripts/pprocess` predates the built-in parallel-processing
support entirely. Instead of using `-P`/`-j`/`-f`, it splits an `.egsinp`
file textwise into `N` numbered copies (`basename_w1.egsinp`,
`basename_w2.egsinp`, ...) with adjusted random-number seeds and history
counts, and submits each as an independent job. It only understands the
classic RZ/BEAM input format (`dosxyznrc`, `BEAM*`) and is only useful today
if you don't have a C compiler available (the built-in parallel support
requires `egs_c_utils.c` to have compiled) or you need to manually split
pieces of a BEAM phase-space source. New work should prefer `-P`/`-j`/`-f`
(via `egs-parallel` or `exb`) instead.

## 4. Combining results and cleaning up

### 4.1 Combining

With the default job-control-file (JCF) run control, the last job to finish
combines all partial results automatically into `basename.egslog` /
`basename.egsdat`. Nothing further is needed in the common case.

If a job died before finishing, or you used a simple RCO (`-s`/`simple`),
combine manually by adding to the `run control` block of the input file and
re-running the application once, non-parallel:

```
:start run control:
    ncase       = ...
    calculation = combine
:stop run control:
```

You can also use `calculation = analyze` to print the stored results from
an existing `.egsdat` file without recomputing anything.

### 4.2 Cleaning up: `egs-parallel-clean`

```
$HEN_HOUSE/scripts/bin/egs-parallel-clean [options] { -f | -n } basename [basename2 ...]
```

| Option | Meaning |
|---|---|
| `-l` / `--list` | just list `.egslog` basenames in the current directory and exit |
| `-f` / `--force` | actually remove temporary files |
| `-n` / `--dry-run` | show what would be removed, without removing it |
| `-x` / `--extra` | also remove the `.egsparallel` log (otherwise it's merged/preserved) |
| `-v` / `--verbose` | log the actual shell commands instead of short messages |

It removes per-run scratch files (`.lock`, `.mederr`, `.eo`, `.e`, `.o`,
`.pbsdsh`, `.egsjob`, work directories `basename_w*`, `egsrun_*_basename_*`)
for each named simulation. `basename` can be a glob pattern in quotes, e.g.
`egs-parallel-clean -f "myinput*"`.

Unless `-x` is given, before deleting it merges the `egs-parallel` log
lines from all involved files into `basename.egsparallel-log`, and all
non-log job output into `basename.egsparallel-out`, so you keep a record of
the run after the scratch files are gone.

### 4.3 `clean_after_parallel` (legacy)

The older, interactive `$HEN_HOUSE/scripts/clean_after_parallel basename`
script is the equivalent cleanup tool for `pprocess`/`exb`-style runs: it
checks that `basename.egslst` is newer than all `basename_w*.egslst` files
(i.e. that the combination step has already run), optionally runs `addphsp`
on split phase-space files, optionally concatenates one worker's `.egslst`
into the main one, and then interactively prompts before deleting the
`basename_w*` files.

## 5. Quick reference

| Task | Modern (`egs-parallel`) | Legacy (`exb`) |
|---|---|---|
| Run N jobs locally | `egs-parallel -b cpu -n N -c '...'` | n/a (use `batch=at` with `now`, or run manually) |
| Run N PBS jobs | `egs-parallel -b pbs -n N -c '...'` | `exb ... batch=pbs p=N` |
| Run as a PBS job array | `egs-parallel -b jobarr -n N -c '...'` | n/a |
| Run via PBS distributed shell | `egs-parallel -b pbsdsh -n N -c '...'` | `exb ... batch=pbsdsh p=N` |
| Run via Slurm (one job per `sbatch`) | `egs-parallel -b slurm -n N -c '...'` | `EGS_BATCH_SYSTEM=slurm exb ... p=N` |
| Run via Slurm (xargs-batched per node, NRC-specific) | n/a | `egs-slurm app inp pegs p=N` (local script, not part of distribution) |
| Resubmit/add jobs to a running set | `egs-parallel ... -f` (force) | `exb ... p=K start=S stop=E fresh=no` |
| Clean up scratch files | `egs-parallel-clean -f basename` | `clean_after_parallel basename` |

See also: PIRS-877 (parallel processing with the NRC user codes) and the
*common inputs* section of PIRS-898 (`run control` block, `-P`/`-j`/`-s`
command-line options) for the application-side details.
