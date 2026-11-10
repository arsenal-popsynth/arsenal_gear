"""
record.py
===========

Store individual timestepped outputs.
"""

from astropy.units import Quantity


def add_record(f, name, time: Quantity["time"], data: Quantity, record_units=None):
    """
    Add a record to the data/ group in the hdf5 file.  Each record
    consists of an array of data values and an array of times.
    """
    ds = f.create_dataset(
        f"data/{name}/qty",
        data.to_value(data.unit if record_units is None else record_units),
    )
    ds.attrs("units", data.unit.name)
    tds = f.create_dataset(f"data/{name}/time", time.value)
    tds.attrs("units", time.unit.name)
