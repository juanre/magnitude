#!/usr/bin/env python

import pytest

import magnitude
from magnitude import Magnitude, mg, new_mag


def test_unit():
    m = mg(1, "1/s")
    assert (m.val, m.unit) == (1.0, [0, -1, 0, 0, 0, 0, 0, 0, 0])

    m = mg(10, "m2 K / min")
    assert m.val - 10.0 / 60.0 < 1e-6
    assert m.unit == [2, -1, 1, 0, 0, 0, 0, 0, 0]
    m = mg(10, "Gib")
    assert m.val == 10737418240.0
    a = mg(1000, "m/s")
    b = a.ounit("km/s")
    assert b.toval() == 1.0


def test_newmag():
    lb = Magnitude(0.45359237, kg=1)
    new_mag("lb", lb)
    assert lb.val - 0.45359237 < 1e-6
    assert lb.unit == [0, 0, 0, 1, 0, 0, 0, 0, 0]
    assert mg(1, "lb") == lb


def test_construct():
    m = mg(10, "m/s") + (1, "lightyear/year")
    n = mg(10, "m/s") + mg(1, "lightyear/year")
    o = mg(1, "lightyear/year") + mg(10, "m/s")
    p = mg(5, "m/s") + mg(1, "lightyear/year") + (5, "m/s")
    assert m == n == o == p and m.unit == [1, -1, 0, 0, 0, 0, 0, 0, 0]

    m = mg(10, "km/h") - (1, "km/h")
    n = mg(10, "km/h") - mg(1, "km/h")
    assert m.unit == [1, -1, 0, 0, 0, 0, 0, 0, 0]
    assert m.val == 10000.0 / 3600 - 1000.0 / 3600
    assert m == n


def test_arithmetic():
    m = mg(10, "m/s") * (2, "s")
    n = mg(20, "m")
    assert m == n
    assert (10, "m/s") * mg(2, "s") == n
    assert mg(10, "m/s") / (0.5, "m") == mg(20, "1/s")
    assert mg(10, "m/s") / mg(0.5, "m") == mg(20, "1/s")
    assert mg(5, "kg/m2") % 3 == mg(2, "kg/m2")
    assert mg(5, "m/s") // mg(3, "m") == mg(1, "1/s")
    assert mg(10, "m/s") ** 2 == mg(100, "m2/s2")
    assert mg(100, "m2/s2") ** 0.5 == mg(10, "m/s")
    assert (-mg(10, "m/s")).val == -10.0
    assert (+mg(10, "m/s")).val == 10.0
    assert abs(-mg(10, "m/s")) == mg(10, "m/s")
    assert int(mg(10.9, "m/s")) == 10
    assert float(mg(10.9, "m/s")) == 10.9
    assert int(mg(10.9, "m/s")) == 10

    assert mg(10, "s-12 / m-12") / mg(0.5, "m") == mg(20, "m11 /s12")
    assert mg(10, "s-12 / m-12 cd2") / mg(0.5, "m") == mg(20, "m11 cd-2/s12")


def test_selfmod():
    m = mg(10, "m/s")
    m += mg(3, "m/s")
    assert m.val == 13.0
    m -= mg(43, "m/s")
    assert m.val == -30.0
    m *= mg(2, "s")
    assert m == mg(-60, "m")

    m = mg(10, "m2/s")
    m /= (5, "m")
    assert m == mg(2, "m/s")

    m = mg(5, "m/s")
    m //= mg(3, "m")
    assert m == mg(1, "1/s")

    m = mg(5, "m/s")
    m = m**2
    assert m == mg(25, "m2/s2")


def test_comp():
    assert mg(1, "m/s") != mg(2, "m/s")
    assert mg(1, "m/s") < mg(2, "m/s")
    assert not mg(1, "m/s") > mg(2, "m/s")


def test_scalar():
    assert (2 * mg(10, "m/s")) == mg(20, "m/s")


def test_coerce_edge_cases():
    """Test edge cases in coerce method."""
    m = mg(1, "m")

    # Test single-element tuple
    result = m.coerce((5.0,))
    assert result.val == 5.0
    assert result.dimensionless()

    # Test invalid tuple length
    with pytest.raises(magnitude.ConversionError) as exc_info:
        m.coerce((1, 2, 3))
    assert "Cannot coerce tuple of length 3" in str(exc_info.value)

    # Test non-numeric type
    with pytest.raises(magnitude.ConversionError) as exc_info:
        m.coerce("string")
    assert "Cannot coerce str to Magnitude" in str(exc_info.value)


def test_resolution_parsing_error():
    """Test error case in resolution parsing."""
    # Test invalid resolution string
    with pytest.raises(magnitude.ConversionError) as exc_info:
        magnitude._parse_resolution_dimensions("invalid_resolution")
    assert "Cannot parse resolution string" in str(exc_info.value)


def test_ensmg_edge_cases():
    """Test edge cases in ensmg function."""
    # Test single-element tuple with unit
    result = magnitude.ensmg((5.0,), "m")
    assert result == mg(5.0, "m")

    # Test single-element tuple without unit
    result = magnitude.ensmg((5.0,))
    assert result.val == 5.0
    assert result.dimensionless()

    # Test invalid tuple length
    with pytest.raises(magnitude.ConversionError) as exc_info:
        magnitude.ensmg((1, 2, 3))
    assert "Can't convert (1, 2, 3) to Magnitude" in str(exc_info.value)

    # Test non-convertible type
    with pytest.raises(magnitude.ConversionError) as exc_info:
        magnitude.ensmg("not a number")
    assert "Can't convert" in str(exc_info.value)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
