"""Test reverse operations for Magnitude."""

import pytest

from magnitude import IncompatibleUnitsError, mg


def test_radd():
    """Test reverse addition."""
    # Scalar + Magnitude
    assert 5 * mg(10, "m") == mg(50, "m")

    # Tuple + Magnitude
    assert (5, "m") + mg(10, "m") == mg(15, "m")

    # Should fail with incompatible units
    with pytest.raises(IncompatibleUnitsError):
        5 + mg(10, "m")


def test_rsub():
    """Test reverse subtraction."""
    # Tuple - Magnitude
    assert (15, "m") - mg(10, "m") == mg(5, "m")

    # Should fail with incompatible units
    with pytest.raises(IncompatibleUnitsError):
        5 - mg(10, "m")

    # Should fail with incompatible units
    with pytest.raises(IncompatibleUnitsError):
        (15, "kg") - mg(10, "m")


def test_rmul():
    """Test reverse multiplication."""
    # Scalar * Magnitude
    assert 5 * mg(10, "m") == mg(50, "m")

    # Tuple * Magnitude
    assert (5, "s") * mg(10, "m") == mg(50, "m s")

    # Complex number * Magnitude
    result = (2 + 3j) * mg(10, "m")
    assert result.val == (20 + 30j)
    assert result.unit == [1, 0, 0, 0, 0, 0, 0, 0, 0]


def test_rtruediv():
    """Test reverse true division."""
    # Scalar / Magnitude
    assert 10 / mg(2, "m") == mg(5, "1/m")

    # Tuple / Magnitude
    assert (10, "m") / mg(2, "s") == mg(5, "m/s")

    # Check units work correctly
    result = 1 / mg(2, "m/s")
    assert result.has_dimension("s/m")


def test_rfloordiv():
    """Test reverse floor division."""
    # Scalar // Magnitude
    assert 25 // mg(10, "m") == mg(2, "1/m")

    # Tuple // Magnitude
    assert (25, "m") // mg(10, "s") == mg(2, "m/s")

    # Test with negative numbers
    assert -25 // mg(10, "m") == mg(-3, "1/m")


def test_rmod():
    """Test reverse modulo."""
    # Scalar % Magnitude
    result = 25 % mg(10, "m")
    assert result.val == 5.0
    assert result.dimensionless()

    # Tuple % Magnitude
    assert (25, "m") % mg(3, "m") == mg(1, "m")

    # Test edge case
    assert (10, "m") % mg(10, "m") == mg(0, "m")


def test_rdivmod():
    """Test reverse divmod."""
    # Scalar divmod
    quot, rem = divmod(25, mg(10, "m"))
    assert quot == mg(2, "1/m")
    assert rem.val == 5.0

    # Tuple divmod
    quot, rem = divmod((25, "m"), mg(10, "m"))
    assert quot.val == 2.0
    assert quot.dimensionless()
    assert rem == mg(5, "m")


def test_edge_cases():
    """Test edge cases for reverse operations."""
    # Zero values
    assert 0 * mg(10, "m") == mg(0, "m")
    assert (0, "m") + mg(10, "m") == mg(10, "m")

    # Negative values
    assert -5 * mg(10, "m") == mg(-50, "m")
    assert (-5, "m") - mg(10, "m") == mg(-15, "m")

    # Float precision
    result = 0.1 * mg(10, "m")
    assert abs(result.val - 1.0) < 1e-10


def test_chained_operations():
    """Test chained reverse operations."""
    # Multiple operations
    result = 2 * mg(5, "m") + (3, "m")
    assert result == mg(13, "m")

    # Complex expression
    result = 100 / (2 * mg(5, "m/s"))
    assert result == mg(10, "s/m")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
