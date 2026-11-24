"""
record.py
===========

Store individual timestepped outputs.
"""

import h5py
from astropy.units import Quantity


def write_recordset(f, records: dict):
    """
    Write a complete set of records into the hdf5 file f, using add_record for
    each one.

    The function does a first pass over all the records to find the non-unique
    time arrays, and then replaces all duplicates with strings indicating the
    first instance of that particular time array.  That then allows add_record
    to write the first instance as an actual array dataset, and then subsequent
    ones as softlinks.  Time arrays only count as duplicates if they have the
    same units, shape, dtype, and values once converted to the units they will
    be written in.

    :param f: Open, writable hdf5 file to write the records into
    :type f: h5py.File
    :param records: Records to write, keyed by record name.  Each value is a
        dict with required keys ``data`` (a Quantity array of values) and
        ``time`` (a Quantity array of times, or the name of another record to
        share times with), and optional keys ``record_units``, ``time_units``
        (astropy units to convert to before writing) and ``descriptor`` (a
        short string describing the record, default ``""``).
    :type records: dict
    :raises TypeError: If any record's time is not a Quantity or a str
    :raises astropy.units.UnitConversionError: If a requested unit is
        incompatible with the record's data or times
    """
    first_seen = {}  # the time values as they would be written -> record name
    times = {}
    for name, record in records.items():
        time = record["time"]
        if isinstance(time, Quantity):
            unit = record.get("time_units")
            unit = time.unit if unit is None else unit
            values = time.to_value(unit)
            # This key is a bit more complicated than seems necessary because
            # the tobytes() return could in rare occasions cause collisions if
            # two arrays have different shape AND different dtypes
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
            descriptor=record.get("descriptor", ""),
        )


def add_record(
    f, name, time, data: Quantity, record_units=None, time_units=None, descriptor=""
):
    """
    Add a record to the data/ group in the hdf5 file.  Each record
    consists of an array of data values and an array of times, written to
    ``data/<name>/qty`` and ``data/<name>/time``.

    The times can either be a Quantity array, which is written as a new
    dataset, or a str naming another record, in which case a SoftLink to that
    record's times is written instead (and time_units is ignored).  The linked
    record does not need to have been written yet.

    :param f: Open, writable hdf5 file to write the record into
    :type f: h5py.File
    :param name: Name of the record, used as its group name under ``data/``
    :type name: str
    :param time: Times the data values apply at, or the name of another
        record whose times should be shared via a SoftLink
    :type time: Quantity["time"] or str
    :param data: Values of the record
    :type data: Quantity
    :param record_units: Units to convert data to before writing (defaults to
        the units of data)
    :type record_units: astropy.units.Unit, optional
    :param time_units: Units to convert time to before writing (defaults to
        the units of time; ignored if time is a str)
    :type time_units: astropy.units.Unit, optional
    :param descriptor: Short human-readable description of the record
    :type descriptor: str, optional
    :raises TypeError: If time is not a Quantity or a str
    :raises astropy.units.UnitConversionError: If record_units or time_units
        is incompatible with data or time
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
    ds.attrs["descriptor"] = descriptor

    if isinstance(time, str):
        f[f"data/{name}/time"] = h5py.SoftLink(f"/data/{time}/time")
    else:
        tds = f.create_dataset(
            f"data/{name}/time",
            data=time.to_value(time.unit if time_units is None else time_units),
        )
        tds.attrs["units"] = time.unit.name if time_units is None else time_units.name
