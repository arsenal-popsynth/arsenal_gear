"""
record.py
===========

Store individual timestepped outputs.
"""

import h5py
from astropy.units import Quantity


def add_record(f, name, time, data: Quantity, record_units=None, time_units=None):
    """
    Add a record to the data/ group in the hdf5 file.  Each record
    consists of an array of data values and an array of times.

    The times can either be a Quantity array, which is written as a new
    dataset, or a str naming another record, in which case a SoftLink to that
    record's times is written instead (and time_units is ignored).
    """
    if not isinstance(time, (Quantity, str)):
        raise TypeError(
            f"time must be a Quantity or a str naming another record, not {type(time)}"
        )

    ds = f.create_dataset(
        f"data/{name}/qty",
        data=data.to_value(data.unit if record_units is None else record_units),
    )
    ds.attrs["units"] = data.unit.name if record_units is None else record_units.name

    if isinstance(time, str):
        f[f"data/{name}/time"] = h5py.SoftLink(f"/data/{time}/time")
    else:
        tds = f.create_dataset(
            f"data/{name}/time",
            data=time.to_value(time.unit if time_units is None else time_units),
        )
        tds.attrs["units"] = time.unit.name if time_units is None else time_units.name
