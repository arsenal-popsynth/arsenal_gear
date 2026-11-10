"""
record.py
===========

Store individual timestepped outputs.
"""

import h5py
from astropy.units import Quantity


def write_recordset(f, records: dict):
    """
    This function takes a dictionary of records and uses add_record to write
    them into the hdf5 file f.  Each record's key is the name of the record,
    and the value for that key is itself a two-item dict.  The 'data' item
    contains the data, and the 'time' unit contains the time, both of which
    will get get passed to add record.  They may also include optional
    record_units and time_units keys.

    The function does a first pass over all the records to find the non-unique
    time arrays, and then replaces all duplicates with strings indicating the
    first instance of that particular time array.  That then allows add_record
    to write the first instance as an actual array dataset, and then subsequent
    ones as softlinks.
    """
    first_seen = {}  # the time values as they would be written -> record name
    times = {}
    for name, record in records.items():
        time = record["time"]
        if isinstance(time, Quantity):
            unit = record.get("time_units")
            unit = time.unit if unit is None else unit
            values = time.to_value(unit)
            key = (unit.name, values.shape, values.dtype.str, values.tobytes())
            if key in first_seen:
                time = first_seen[key]
            else:
                first_seen[key] = name
        times[name] = time

    for name, record in records.items():
        add_record(
            f,
            name,
            times[name],
            record["data"],
            record_units=record.get("record_units"),
            time_units=record.get("time_units"),
        )


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
