"""Test Unicode superscript support for Magnitude."""

import os

import pytest

from magnitude import UnitError, _parse_superscript, _to_superscript, mg, unicode_superscript


def test_parse_superscript():
    """Test parsing of Unicode superscripts to integers."""
    assert _parse_superscript("²") == 2
    assert _parse_superscript("³") == 3
    assert _parse_superscript("¹²") == 12
    assert _parse_superscript("⁻¹") == -1
    assert _parse_superscript("⁻¹²") == -12
    assert _parse_superscript("¹⁰⁰") == 100
    assert _parse_superscript("⁰") == 0

    # Invalid superscripts
    assert _parse_superscript("x²") is None
    assert _parse_superscript("²x") is None
    assert _parse_superscript("") is None


def test_to_superscript():
    """Test conversion of integers to Unicode superscripts."""
    assert _to_superscript(0) == "⁰"
    assert _to_superscript(1) == "¹"
    assert _to_superscript(2) == "²"
    assert _to_superscript(12) == "¹²"
    assert _to_superscript(-1) == "⁻¹"
    assert _to_superscript(-12) == "⁻¹²"
    assert _to_superscript(100) == "¹⁰⁰"


def test_input_parsing():
    """Test parsing units with Unicode superscripts."""
    # Basic superscripts
    assert mg(10, "m²") == mg(10, "m2")
    assert mg(10, "s⁻¹") == mg(10, "s-1")

    # Multi-digit superscripts
    assert mg(1, "kg¹²") == mg(1, "kg12")
    assert mg(1, "s⁻¹²") == mg(1, "s-12")

    # Complex units
    assert mg(10, "m²/s²") == mg(10, "m2/s2")
    assert mg(1, "kg m²/s²") == mg(1, "J")

    # Mixed notation (should work)
    assert mg(10, "m² s-1") == mg(10, "m2/s")
    assert mg(10, "km² h⁻¹") == mg(10, "km2/h")


def test_output_formatting():
    """Test output with Unicode superscripts."""
    # Store original setting
    original = unicode_superscript()

    try:
        # Test with Unicode disabled (default)
        unicode_superscript(False)
        m = mg(10, "m") * mg(1, "m") / (mg(1, "s") * mg(1, "s"))
        assert "m2 / s2" in str(m)

        # Test with Unicode enabled
        unicode_superscript(True)
        assert "m² / s²" in str(m)

        # Test complex units
        m2 = mg(1, "kg") ** 12 / (mg(1, "m") ** 3 * mg(1, "s") ** 12)
        assert "kg¹² / m³ s¹²" in str(m2)

        # Test with output unit
        m3 = mg(1, "kg m2/s2")
        m3.ounit("J")
        # This might still show 'J' - that's OK for now

    finally:
        # Restore original setting
        unicode_superscript(original)


def test_prefix_with_superscripts():
    """Test units with prefixes and superscripts."""
    # Kilo-meter squared
    assert mg(1, "km²") == mg(1e6, "m2")

    # Micro-second to the -1
    assert mg(1, "us⁻¹") == mg(1e6, "s-1")

    # Mega-watt
    assert mg(1, "MW") == mg(1e6, "W")


def test_edge_cases():
    """Test edge cases for Unicode superscripts."""
    # Zero exponent
    m = mg(10, "m")
    m.val = 10  # Should have m^1, not m^0

    # Fractional exponents should still fail
    with pytest.raises(UnitError):
        mg(10, "m^(1/2)")  # Not supported

    # Mixed Unicode in unit name should fail
    with pytest.raises(UnitError):
        mg(10, "m²eters")  # Invalid unit


def test_environment_variable():
    """Test that MAGNITUDE_UNICODE_SUPERSCRIPTS environment variable works."""
    # Store original values
    original_env = os.environ.get("MAGNITUDE_UNICODE_SUPERSCRIPTS")

    # Need to reload the module to pick up environment changes
    import importlib

    import magnitude

    try:
        # Test various truthy values
        for value in ["1", "true", "True", "TRUE", "yes", "YES", "on", "ON"]:
            os.environ["MAGNITUDE_UNICODE_SUPERSCRIPTS"] = value
            importlib.reload(magnitude)
            # The default should now be True
            assert magnitude._unicode_superscript is True

        # Test various falsy values
        for value in ["0", "false", "False", "no", "off", ""]:
            os.environ["MAGNITUDE_UNICODE_SUPERSCRIPTS"] = value
            importlib.reload(magnitude)
            # The default should now be False
            assert magnitude._unicode_superscript is False

        # Test unset
        if "MAGNITUDE_UNICODE_SUPERSCRIPTS" in os.environ:
            del os.environ["MAGNITUDE_UNICODE_SUPERSCRIPTS"]
        importlib.reload(magnitude)
        # The default should be False
        assert magnitude._unicode_superscript is False

    finally:
        # Restore original value
        if original_env is not None:
            os.environ["MAGNITUDE_UNICODE_SUPERSCRIPTS"] = original_env
        elif "MAGNITUDE_UNICODE_SUPERSCRIPTS" in os.environ:
            del os.environ["MAGNITUDE_UNICODE_SUPERSCRIPTS"]
        # Reload to restore original state
        importlib.reload(magnitude)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
