# Arsenal Output Formats
## HDF5 Outputs

Arsenal Gear writes its binary outputs as [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
files using [h5py](https://docs.h5py.org/). Each file has two top-level groups:

- `metadata/` records how and where the file was made, so you can reproduce the
  run from the file alone.
- `data/` holds the time series data ("records") the run produced.  These can be
  either population-averaged or per star for discrete stellar populations.

```text
output.h5
├── metadata/
│   ├── rundate        ISO 8601 timestamp of when the file was written
│   ├── hostname       machine the run was executed on
│   ├── paramfile      full text of the parameter file used for the run
│   ├── bibtex         (optional) BibTeX for the models used
│   └── source         gzipped tarball of the arsenal_gear source tree
│                        attrs: format, git-commit, git-dirty
└── data/
    ├── <record name>/
    │   ├── qty        the record's values           attrs: units, descriptor
    │   └── time       the times those values apply at   attrs: units
    └── ...
```

### Metadata

{py:func}`arsenal_gear.file_io.metadata.construct_meta_dataset` fills in the `metadata/`
group:

| Dataset              | Contents                                                      |
| -------------------- | ------------------------------------------------------------- |
| `metadata/rundate`   | `datetime.now().isoformat()` at write time                    |
| `metadata/hostname`  | `socket.gethostname()`                                        |
| `metadata/paramfile` | The parameter file text, passed in as a string                |
| `metadata/bibtex`    | BibTeX string for the run, written only if one is given       |
| `metadata/source`    | The installed `arsenal_gear` package as a `.tar.gz` (opaque bytes) |

The `source` dataset carries three attributes:

- `format`: always `"tar.gz"`.
- `git-commit`: the commit hash of the source tree, or `""` if the package isn't a git
  checkout or `git` isn't installed.
- `git-dirty`: `True` if the working tree had uncommitted changes. When no commit can be
  found this is also `True`, since the source can't be matched to a known commit.

This source tarball is __an exact copy of the source that generated the hdf5
file__.  It doesn't matter how or where you have installed arsenal gear, that
source tarball contains the source that made the hdf5 file: re-running with that
copy of the source will produce a bitwise-identical result for everything except
the `metadata/` (as `hostname` and `rundate` may be different).

The tarball leaves out every directory named `data` or `__pycache__`, which keeps
large downloaded model grids out of each output file. The source of the run is kept
in the file, but the stellar evolution data it read is not.

To recover the source from an output file, use the `arsenal-extract-source` command
that is installed with Arsenal Gear:

```text
$ arsenal-extract-source output.h5 --info     # print provenance only
   rundate: 2026-11-24T16:52:31.293909
  hostname: mymachine
git-commit: 6117e5201f83f89c9563242fc4bd940d5065f9d0
 git-dirty: False
$ arsenal-extract-source output.h5            # extract to ./output_source/
$ arsenal-extract-source output.h5 -o src/    # extract to ./src/
```

It won't extract into a directory that already contains files unless you pass `--force`.

### Records

A *record* is a single time series, such as the supernova energy injected as a function of
population age. Each record is stored as a group under `data/` holding two datasets:

- `qty`: the values, with a `units` attribute (an astropy unit name such as `"erg"`) and
  a `descriptor` attribute (a short human-readable description).
- `time`: the times, with its own `units` attribute.

Values are stored as plain numbers. The units travel alongside them as attributes, so
you can rebuild an astropy {py:class}`~astropy.units.Quantity` when reading:

```python
import astropy.units as u
import h5py

with h5py.File("output.h5", "r") as f:
    qty = f["data/sn_energy/qty"]
    time = f["data/sn_energy/time"]
    energy = qty[()] * u.Unit(qty.attrs["units"])
    age = time[()] * u.Unit(time.attrs["units"])
    print(qty.attrs["descriptor"])
```

#### Shared time arrays

Many records are sampled at the same times. Rather than storing an identical `time`
array for each one, Arsenal Gear writes the array once and makes the other records'
`time` entries HDF5 soft links to it. Reading through h5py follows the link for you, so
`f["data/<name>/time"]` works the same either way. To check whether a record's times are
a link:

```python
link = f["data/wind_mass"].get("time", getlink=True)
if isinstance(link, h5py.SoftLink):
    print("times shared with", link.path)  # e.g. /data/sn_energy/time
```

Two time arrays count as duplicates only if they match exactly once converted to the
units they will be written in: same unit, shape, dtype, and values. Arrays that are
physically equal but stored in different units (say `Myr` and `yr`) are written
separately.

### Writing records

{py:func}`arsenal_gear.file_io.record.write_recordset` is the main entry point. It takes a
dict that maps each record name to a dict with these keys:

| Key            | Required | Meaning                                                   |
| -------------- | -------- | --------------------------------------------------------- |
| `data`         | yes      | `Quantity` array of values                                |
| `time`         | yes      | `Quantity` array of times, or the name of another record to share times with |
| `record_units` | no       | Unit to convert `data` to before writing (defaults to `data.unit`) |
| `time_units`   | no       | Unit to convert `time` to before writing (defaults to `time.unit`) |
| `descriptor`   | no       | Short description stored on `qty` (defaults to `""`)      |

```python
import astropy.units as u
import h5py
import numpy as np

from arsenal_gear.file_io.metadata import construct_meta_dataset
from arsenal_gear.file_io.record import write_recordset

age = np.linspace(0, 40, 100) * u.Myr

with h5py.File("output.h5", "w") as f:
    construct_meta_dataset(f, paramfile=open("my_params.toml").read())
    write_recordset(
        f,
        {
            "sn_energy": {
                "data": e_sn,  # Quantity in erg
                "time": age,
                "descriptor": "Cumulative SN energy",
            },
            "wind_mass": {
                "data": m_wind,  # Quantity in Msun
                "time": age,  # identical to sn_energy, so it becomes a soft link
                "record_units": u.g,
                "descriptor": "Cumulative wind mass loss",
            },
        },
    )
```

`write_recordset` first scans every record to find duplicate time arrays, then writes
each record with {py:func}`~arsenal_gear.file_io.record.add_record`. Only the first record
with a given time array stores it; later duplicates get a soft link to it.

You can also call `add_record` directly to write a single record. Pass a record name as
`time` to link to that record's times explicitly. The link target doesn't have to exist
yet; it resolves once that record is written.

```python
from arsenal_gear.file_io.record import add_record

add_record(f, "sn_energy", age, e_sn, descriptor="Cumulative SN energy")
add_record(f, "sn_mass", "sn_energy", m_sn, descriptor="Cumulative SN ejecta mass")
```

If you link this way, `time_units` is ignored and the record inherits the linked times'
units. Any other `time` type raises a `TypeError`, and incompatible units raise
{py:class}`astropy.units.UnitConversionError`.

## API Reference

```{eval-rst}
.. autofunction:: arsenal_gear.file_io.metadata.construct_meta_dataset
   :no-index:

.. autofunction:: arsenal_gear.file_io.metadata.get_full_src_tree
   :no-index:

.. autofunction:: arsenal_gear.file_io.record.write_recordset
   :no-index:

.. autofunction:: arsenal_gear.file_io.record.add_record
   :no-index:
```
