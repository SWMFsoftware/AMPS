# Qualification-machine profiles

Every invocation of `src/test/run_tests.py` must explicitly select at least
one machine profile with `--config`. This prevents a long or distributed run
from silently using workstation defaults, an unintended MPI implementation,
or repository-local storage.

Profiles are named after `hostname -s`:

- `passion.json` is the maintained profile for the Passion host;
- `template.json` is a copy-and-edit starting point for another machine.

Create a profile for a new host with:

```sh
cp test/qualification-machine/template.json \
  test/qualification-machine/$(hostname -s).json
```

Edit every path and tool command before running tests. In particular, choose a
durable `output_root`, a node-local `temporary_root`, the reference-data root,
SPICE installation and kernels, compiler/MPI wrappers, rank spelling, and
resource limits. Remove optional tool entries that are not installed; the
corresponding portability identities will then report `SKIP` rather than
attempting an unrelated executable from `PATH`.

Run on Passion with:

```sh
python3 src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --campaign nightly --streams 3 --jobs 2
```

`--streams` bounds the independent core, Enceladus, and Moon execution
streams. `--jobs` separately bounds nonexclusive suites inside each stream.
Native/configuration-changing suites remain globally serialized regardless of
either value.

The make wrapper chooses the matching short-hostname file itself:

```sh
make -C src/test nightly
make -C src/test CONFIG=/path/to/reviewed-machine.json nightly
make -C src/test STREAMS=3 JOBS=2 nightly
```

Direct Python invocations never choose a profile implicitly; this make
behavior works because the wrapper visibly adds its own `--config` argument.
Before launching a costly campaign, validate a new profile and inspect the
resolved registry without running tests:

```sh
python3 -m json.tool test/qualification-machine/$(hostname -s).json
python3 src/test/run_tests.py \
  --config test/qualification-machine/$(hostname -s).json \
  --list --streams 3
```

Run the complete expensive/qualification selection with:

```sh
python3 src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --full-validation --streams 3 --jobs 2
```

The committed default file remains the schema and fallback-value base; a
qualification-machine file is an explicit overlay. Multiple `--config` options
are allowed and are applied from left to right, but at least one is mandatory.
Do not place credentials in these files.
