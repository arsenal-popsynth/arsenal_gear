"""
element_yields
==========

Loaders and functions for element yields.
"""

from .limongichieffi2018 import LimongiChieffi2018
from .nugrid import Battino20192021, NuGrid, Pignatari2016, Ritter2018
from .sukhbold2016 import Sukhbold2016

# fmt: off
elements = [ "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
"Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
"Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
"Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",
"I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
"Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
"Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
"U"]
# fmt: on

__all__ = [
    "LimongiChieffi2018",
    "Battino20192021",
    "NuGrid",
    "Pignatari2016",
    "Ritter2018",
    "Sukhbold2016",
]
