"""
record.py
===========

Store individual timestepped outputs.
"""

from astropy.units import Quantity


def add_record(f, name, time: Quantity["time"], data: Quantity):
    """
    Add a record to the data/ group in the hdf5 file.  Each record
    consists of an array of data values and an array of times.
    """
    ds = f.create_dataset(f"data/{name}/qty", data.value)
    ds.attrs("units", data.unit.name)
    tds = f.create_dataset(f"data/{name}/time", time.value)
    tds.attrs("units", time.unit.name)
