import pytest

import arsenal_gear

def test_star():
    star = arsenal_gear.star.Star(mass=1.0,metals=0.02,age=1.0)
    assert star.mass == 1.0
    assert star.metals == 0.02
    assert star.age == 1.0

def test_binary():
    star1 = arsenal_gear.star.Star(mass=1.0,metals=0.02,age=1.0)
    star2 = arsenal_gear.star.Star(mass=2.0,metals=0.01,age=0.0)
    binary = arsenal_gear.binary.Binary(star_1=star1, star_2=star2, radius=10, eccentricity=0.1)
    # Check the first star
    assert binary.star_1.mass == star1.mass
    assert binary.star_1.metals == star1.metals
    assert binary.star_1.age == star1.age

    # Check the second star
    assert binary.star_2.mass == star2.mass
    assert binary.star_2.metals == star2.metals
    assert binary.star_2.age == star2.age

    # Check binary properties
    assert binary.radius == 10
    assert binary.eccentricity == 0.1
