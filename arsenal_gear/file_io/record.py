"""
record.py
===========

Store individual timestepped outputs.
"""

import h5py
from astropy.units import Quantity


class TimeRecord:
    """
    The times associated with a record can be in one of three possible forms:
    an array of times (the simplest case); a set of limits plus a step size in
    log or linear space; or a SoftLink to another time record.
    """

    def __init__(self, time):
        if isinstance(time, Quantity):
            self.type = "array"
            self.time = time
        elif isinstance(time, str):
            self.type = "link"
            self.time = time
        else:
            raise TypeError(
                f"time must be a Quantity or a str naming another record, not {type(time)}"
            )

    def write(self, f, name: str, time_units=None):
        if self.type == "array":
            self._array_write(f, name, time_units)
        if self.type == "link":
            self._link_write(f, name)

    def _link_write(self, f, name):
        f[f"data/{name}/time"] = h5py.SoftLink(f"/data/{self.time}/time")

    def _array_write(self, f, name: str, time_units=None):
        tds = f.create_dataset(
            f"data/{name}/time",
            data=self.time.to_value(
                self.time.unit if time_units is None else time_units
            ),
        )
        tds.attrs["units"] = (
            self.time.unit.name if time_units is None else time_units.name
        )


def add_record(
    f, name, time: TimeRecord, data: Quantity, record_units=None, time_units=None
):
    """
    Add a record to the data/ group in the hdf5 file.  Each record
    consists of an array of data values and an array of times.
    """
    ds = f.create_dataset(
        f"data/{name}/qty",
        data=data.to_value(data.unit if record_units is None else record_units),
    )
    ds.attrs["units"] = data.unit.name if record_units is None else record_units.name
    time.write(f, name, time_units)
