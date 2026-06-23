# Getting started with parallel EGSnrc jobs

*A quick introduction — see `PARALLEL_JOBS.md` for full details.*

## 1. Why bother?

- EGSnrc simulations are "embarrassingly parallel": every history is
  statistically independent.
- Splitting `N` histories into `N` jobs and combining the results at the
  end gives a near-linear speedup, with no change to your physics setup.
- The tool that does this is **`egs-parallel`**
  (`$HEN_HOUSE/scripts/bin/egs-parallel`), already on your `PATH` if you've
  sourced the standard EGSnrc shell setup.

## 2. Step one: run several jobs on your own computer

No cluster, no scheduler — just use the cores you already have:

```sh
cd $EGS_HOME/egs_chamber
egs-parallel --batch cpu -n 12 -v --command 'egs_chamber -i slab -p 521icru'
```

- `--batch cpu` runs the jobs as local background processes.
- `-n 12` launches 12 copies; pick a number around your core count.
- `--command '...'` is your normal application invocation — `egs-parallel`
  adds the parallel-job flags (`-P`/`-j`/`-f`) for you.
- `-v` prints progress messages to the terminal as it runs.

When all jobs finish, the partial results are combined automatically into
the usual `slab.egslog` / `.egsdat` output — nothing extra to do.

This is the right way to get comfortable with the tool before moving to a
cluster: same command, same output, just `--batch cpu` instead of a
scheduler name.

## 3. Step two: moving to a cluster

Once one machine isn't enough, change `--batch` to match whatever
scheduler your cluster runs — the rest of the command stays the same:

| `--batch` | What it does |
|---|---|
| `pbs` | One PBS job per parallel job (`qsub` × N) |
| `jobarr` | One PBS **job array** covering all N jobs |
| `pdsh` | One PBS job reserving N nodes, fanned out with `pdsh`/`ssh` |
| `pbsdsh` | One PBS job reserving N cores, fanned out with `pbsdsh` |
| `slurm` | One Slurm job per parallel job (`sbatch` × N) |

Example (Slurm):

```sh
egs-parallel --batch slurm -q standard -n 50 \
    -o '--account=myproject --qos=low --time=02:00:00' \
    --command 'egs_chamber -i slab -p 521icru'
```

There's also a legacy front-end, `exb`/`run_user_code_batch`, predating
`egs-parallel`, still installed and usable on sites that already have it
configured — worth knowing it exists, but not the recommended starting
point for new work.

## 4. Where to go next

`PARALLEL_JOBS.md` (same directory) covers, in depth:

- all `egs-parallel` options and batch backends, with more examples;
- the legacy `exb` system and its batch-system definitions;
- how parallel jobs coordinate under the hood (lock file vs. "simple" vs.
  "uniform" run-control objects), and when you'd want to switch;
- combining results manually and cleaning up after a run.
