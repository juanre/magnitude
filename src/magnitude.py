# magnitude  -- a module for computing with numbers with units.
#
# Copyright (c) 2006-2025 Juan Reyero (http://juanreyero.com).
#
# Licensed under the MIT License - see LICENSE file for details.
#
# Home page: http://juanreyero.com/open/magnitude/

"""Python library for computing with physical quantities with units."""

import math
import numbers
import re
from typing import Any, Optional, Union

__version__ = "1.1.0"


class MagnitudeError(Exception):
    """Base exception for magnitude errors."""

    pass


class UnitError(MagnitudeError):
    """Raised when there's an issue with unit parsing or recognition."""

    pass


class IncompatibleUnitsError(MagnitudeError):
    """Raised when operations are attempted on incompatible units."""

    def __init__(self, unit1: list[float], unit2: list[float], operation: str = "operation"):
        self.unit1 = unit1
        self.unit2 = unit2
        self.operation = operation
        super().__init__(f"Incompatible units for {operation}: {unit1} and {unit2}")


class ConversionError(MagnitudeError):
    """Raised when a value cannot be converted to a Magnitude."""

    pass


# Base magnitude names and prefixes.  The _mags dictionary, initialized
# at the end, will contain all the known magnitudes.  Units are
# 9-element arrays, each element the exponent of the unit named by the
# Uname in the same position.

_mags: dict[str, "Magnitude"] = {}
_unames: list[str] = ["m", "s", "K", "kg", "A", "mol", "cd", "$", "b"]
_prefix: dict[str, float] = {
    "y": 1e-24,  # yocto
    "z": 1e-21,  # zepto
    "a": 1e-18,  # atto
    "f": 1e-15,  # femto
    "p": 1e-12,  # pico
    "n": 1e-9,  # nano
    "u": 1e-6,  # micro
    "m": 1e-3,  # mili
    "c": 1e-2,  # centi
    "d": 1e-1,  # deci
    "k": 1e3,  # kilo
    "M": 1e6,  # mega
    "G": 1e9,  # giga
    "T": 1e12,  # tera
    "P": 1e15,  # peta
    "E": 1e18,  # exa
    "Z": 1e21,  # zetta
    "Y": 1e24,  # yotta
    # Binary prefixes, approved by the International
    # Electrotechnical Comission in 1998.  Since then, kb means
    # 1000 bytes; for 1024 bytes use Kib (note the capital K in
    # the binary version, and the lower case for the b of byte,
    # see comment in byte definition below).
    "Ki": 2**10,  # Kibi (<- kilo, 10^3)
    "Mi": 2**20,  # Mebi (<- mega, 10^6)
    "Gi": 2**30,  # Gibi (<- giga, 10^9)
    "Ti": 2**40,  # Tebi (<- tera, 10^12)
    "Pi": 2**50,  # Pebi (<- peta, 10^15)
    "Ei": 2**60,  # Exbi (<- exa, 10^18)
}


# Default print formatting options

_default_prn_format = "%.*f"
_prn_format = _default_prn_format
_prn_prec = 4
_prn_units = True


def reset_default_format() -> None:
    """Resets the default output format.

    By default the output format is "%.*f", where * gets replaced by
    the output precision.
    """
    global _prn_format
    _prn_format = _default_prn_format


def default_format(fmt: Optional[str] = None) -> str:
    """Get or set the default ouptut format.

    Include a fmt if and where you need to specify the output
    precision.  Defaults to %.*f, where the * stands for the
    precision.  Do nothing if fmt is None.

    Returns: default format.

    >>> print(mg(2, 'm2').sqrt())
    1.4142 m
    >>> default_format("%.2f")
    '%.2f'
    >>> print(mg(2, 'm2').sqrt())
    1.41 m
    >>> reset_default_format()
    """
    global _prn_format
    if fmt is not None:
        _prn_format = fmt
    return _prn_format


def output_precision(prec: Optional[int] = None) -> int:
    """Get or set the output precision.

    Package default is 4.  Do nothing is prec is None.

    Returns: default precision.

    >>> default_format("%.*f")
    '%.*f'
    >>> print(mg(2, 'm2').sqrt())
    1.4142 m
    >>> output_precision(6)
    6
    >>> print(mg(2, 'm2').sqrt())
    1.414214 m
    >>> output_precision(4)
    4
    """
    global _prn_prec
    if prec is not None:
        _prn_prec = prec
    return _prn_prec


def output_units(un: Optional[bool] = None) -> bool:
    """Enable or disable the output of units when printing.

    By default output of units is enabled.  Do nothing if un is None.
    When disabled (un is False) print of Magnitudes will produce only
    numbers.

    Return: True if output of units enabled, False otherwise.

    >>> print(mg(2, 'day'))
    2.0000 day
    >>> output_units(False)
    False
    >>> print(mg(2, 'day').ounit('s'))
    172800.0000
    """
    global _prn_units
    if un is not None:
        _prn_units = un
    return _prn_units


def _parse_resolution_dimensions(resolution_str: str) -> tuple[int, int]:
    """Parse a resolution string to extract horizontal and vertical dimensions.

    Args:
        resolution_str: A string in format "[HxV]" or "[N]" where H, V, N are integers.
                       "[HxV]" represents H horizontal x V vertical pixels per inch.
                       "[N]" is shorthand for "[NxN]" (square pixels).

    Returns:
        A tuple of (horizontal_dpi, vertical_dpi)

    Raises:
        ConversionError: If the string cannot be parsed as a resolution.

    Examples:
        >>> _parse_resolution_dimensions("[600x1200]")
        (600, 1200)
        >>> _parse_resolution_dimensions("[300]")
        (300, 300)
    """
    # Try to match HxV pattern
    match = re.search(r"(\d+)x(\d+)", resolution_str)
    if match:
        return int(match.group(1)), int(match.group(2))

    # Check for [N] shorthand format
    if (resolution_str[0] == "[") and (resolution_str[-1] == "]"):
        dpi = int(resolution_str[1:-1])
        return (dpi, dpi)

    raise ConversionError(f"Cannot parse resolution string '{resolution_str}'")


def _is_resolution_string(s: str) -> bool:
    """Check if a string represents a bracketed resolution notation.

    Resolution strings are bracketed values like "[600x600]" or "[300]" used
    in the printing industry to denote pixels per inch.

    Args:
        s: String to check

    Returns:
        True if the string is in bracketed resolution format, False otherwise.
    """
    return (len(s) > 2) and (s[0] == "[") and (s[-1] == "]")


def _resolution_to_square_meters(resolution_str: str) -> float:
    """Convert a resolution string to the area of one pixel in square meters.

    Bracketed resolutions are used in the printing industry to denote
    pixels per inch (DPI). This function calculates the area of a single
    pixel based on the resolution.

    The format can be:
    - "[HxV]" where H and V are horizontal and vertical DPI (e.g., "[300x1200]")
    - "[N]" as shorthand for "[NxN]" for square pixels (e.g., "[600]" = "[600x600]")

    The square brackets indicate we're dealing with area measurements.

    Args:
        resolution_str: Resolution string in bracketed format

    Returns:
        Area of one pixel in square meters

    Examples:
        >>> _resolution_to_square_meters("[600x600]")  # 600 DPI square pixels
        1.792111111111111e-09
        >>> _resolution_to_square_meters("[600]")      # Same as [600x600]
        1.792111111111111e-09
        >>> _resolution_to_square_meters("[150x300]")  # 150x300 DPI pixels
        1.4336888888888889e-08
    """
    h_dpi, v_dpi = _parse_resolution_dimensions(resolution_str)
    # Convert from pixels per inch to meters per pixel
    # 1 inch = 0.0254 meters
    # Area = (0.0254 / h_dpi) * (0.0254 / v_dpi) square meters
    return 0.0254 * 0.0254 / (h_dpi * v_dpi)


def _is_pitch_string(s: str) -> bool:
    """Check if a string represents a pitch notation.

    Pitch strings are dash-delimited values like "-600-" used in the
    printing industry to denote pitch (spacing between elements).

    Args:
        s: String to check

    Returns:
        True if the string is in pitch format (-N-), False otherwise.
    """
    return (len(s) > 2) and (s[0] == "-") and (s[-1] == "-")


def _pitch_to_meters(pitch_str: str) -> float:
    """Convert a pitch string to spacing distance in meters.

    Pitch notation is used in the printing industry to represent spacing
    or pitch between elements. The format "-N-" means 1/N of an inch.

    For example:
    - "-600-" means 1/600 inch pitch (common for dot matrix printers)
    - "-1200-" means 1/1200 inch pitch (finer spacing)

    Args:
        pitch_str: Pitch string in format "-N-" where N is an integer

    Returns:
        The pitch distance in meters

    Examples:
        >>> _pitch_to_meters("-600-")   # 1/600 inch = ~42.3 micrometers
        4.233333333333333e-05
        >>> _pitch_to_meters("-1200-")  # 1/1200 inch = ~21.2 micrometers
        2.1166666666666665e-05
    """
    pitch_value = int(pitch_str[1:-1])  # Extract number between dashes
    # Convert from 1/pitch_value inches to meters
    # 1 inch = 0.0254 meters
    return 0.0254 / pitch_value


class Magnitude:
    """A number with associated units.

    Note: Currently, unit conversions use float arithmetic internally,
    so numeric types like Decimal may lose precision when multiplied
    by units. Future versions may address this limitation.
    """

    def __init__(
        self,
        val: Any,
        m: float = 0,
        s: float = 0,
        K: float = 0,
        kg: float = 0,
        A: float = 0,
        mol: float = 0,
        cd: float = 0,
        dollar: float = 0,
        b: float = 0,
    ) -> None:
        self.val = val
        self.unit: list[float] = [m, s, K, kg, A, mol, cd, dollar, b]
        self.out_unit: Optional[str] = None
        self.out_factor: Optional[Magnitude] = None
        self.oprec: Optional[int] = None
        self.oformat: Optional[str] = None

    def copy(self, with_format: bool = False) -> "Magnitude":
        """Builds and returns a copy of a magnitude.

        The copy includes value and units.  If with_format is set to
        True the default output unit, output factor, output precision
        and output format are also copied.

        >>> a = mg(1000/3., 'mm')
        >>> print(a.output_prec(2))
        333.33 mm
        >>> print(a.copy())
        0.3333 m
        >>> print(a.copy(with_format=True))
        333.33 mm
        """
        cp = Magnitude(self.val, *self.unit)
        if with_format:
            cp.out_unit = self.out_unit
            cp.out_factor = self.out_factor
            cp.oprec = self.oprec
            cp.oformat = self.oformat
        return cp

    def toval(self, ounit: Optional[str] = "") -> Any:
        """Returns the numeric value of a magnitude.

        The value is given in ounit or in the Magnitude's default
        output unit.

        >>> v = mg(100, 'km/h')
        >>> v.toval()
        100.0
        >>> v.toval(ounit='m/s')
        27.77777777777778
        """
        m = self.copy()
        if not ounit:
            ounit = self.out_unit
        if ounit:
            out_factor = self.sunit2mag(ounit)
            m._div_by(out_factor)
        return m.val

    def __repr__(self) -> str:
        """Return a string representation of the Magnitude."""
        return self.__str__()

    def __str__(self) -> str:
        oformat = self.oformat
        oprec = self.oprec
        if oprec is None:
            oprec = _prn_prec
        if oformat is None:
            oformat = _prn_format
        if self.out_unit:
            m = self.copy()
            m._div_by(self.out_factor)
            if "*" in oformat:  # requires the precision arg
                st = oformat % (oprec, m.val)
            else:
                st = oformat % (m.val)
            if _prn_units:
                return st + " " + self.out_unit.strip()
            return st

        if "*" in oformat:
            st = oformat % (oprec, self.val)
        else:
            st = oformat % (self.val)

        if not _prn_units:
            return st

        u = self.unit
        num = " "  # numerator
        for i in range(len(_unames)):
            if u[i] == 1:
                num = num + _unames[i] + " "
            elif u[i] > 0:
                num = num + _unames[i] + str(u[i]) + " "
        den = ""  # denominator
        for i in range(len(_unames)):
            if u[i] == -1:
                den = den + _unames[i] + " "
            elif u[i] < 0:
                den = den + _unames[i] + str(-u[i]) + " "
        if den:
            if num == " ":
                num += "1 "
            st += num + "/ " + den
        elif num != " ":
            st += num
        return st.strip()

    def term2mag(self, s: str) -> "Magnitude":
        """Converts a string with units to a Magnitude.

        Can't divide: use with the numerator and the denominator
        separately (hence the "term").  Returns the Magnitude that the
        string represents.  Units are separated by spaces, powers are
        integers following the unit name.

        Cannot parse fractional units.  Cannot parse multi-digit
        exponents.

        >>> a = mg(1, '')
        >>> print(a.term2mag('V2  A'))
        1.0000 m4 kg2 / s6 A
        >>> print(a.term2mag('kft year')) # kilo-feet year
        9618551037.0820 m s
        """
        m = Magnitude(1.0)
        units = re.split(r"\s", s)
        for u in units:
            if re.search(r"[^\s]", u):
                exp = 1
                if re.search(r"\d$", u):
                    exp = int(u[-1])
                    u = u[0:-1]
                if u in _mags:
                    u = _mags[u].copy()
                elif (len(u) >= 3) and u[0:2] in _prefix and u[2:] in _mags:
                    pr = _prefix[u[0:2]]
                    u = _mags[u[2:]].copy()
                    u.val = pr * u.val
                elif (len(u) >= 2) and u[0] in _prefix and u[1:] in _mags:
                    pr = _prefix[u[0]]
                    u = _mags[u[1:]].copy()
                    u.val = pr * u.val
                elif _is_resolution_string(u):
                    u = Magnitude(_resolution_to_square_meters(u), m=2)
                elif _is_pitch_string(u):
                    u = Magnitude(_pitch_to_meters(u), m=1)
                elif u == "":
                    u = Magnitude(1.0)
                else:
                    raise UnitError(f"Don't know about unit '{u}'")
                for _ in range(exp):
                    m._mult_by(u)
        return m

    def sunit2mag(self, unit: str = "") -> "Magnitude":
        """Convert a units string to a Magnitude.

        Uses term2mag to convert a string with units, possibly
        including a / to separate a numerator and a denominator, to a
        Magnitude.

        >>> a = mg(1, '')
        >>> a.sunit2mag('m/s').toval()
        1.0
        >>> a.sunit2mag('km/h').toval()
        0.2777777777777778
        >>> print(a.sunit2mag('W h'))
        3600.0000 m2 kg / s2
        >>> print(a.sunit2mag('W h').ounit('J'))
        3600.0000 J
        >>> print(a.sunit2mag('m2 kg / s3 Pa'))
        1.0000 m3 / s
        >>> print(a.sunit2mag('m2 kg/s3').ounit('W'))
        1.0000 W
        """
        m = Magnitude(1.0)
        if unit:
            q = re.split(r"/", unit)
            if re.search(r"[^\s]", q[0]):
                m._mult_by(self.term2mag(q[0]))
            if (len(q) == 2) and re.search(r"[^\s]", q[1]):
                m._div_by(self.term2mag(q[1]))
        return m

    def dimensionless(self) -> bool:
        """True if the magnitude's dimension exponents are all zero.

        >>> mg(2, 'K').dimensionless()
        False
        >>> mg(2, 'rad').dimensionless()
        True
        """
        return self.unit == [0] * 9

    def dimension(self) -> list[float]:
        """Return the dimension of the unit in internal (array) format.

        >>> mg(2, 'J').dimension()
        [2, -2, 0, 1, 0, 0, 0, 0, 0]
        """
        return self.unit[:]

    def has_dimension(self, u: str) -> bool:
        """Returns true if the dimension of the magnitude matches u:

        >>> s = mg(120, 'km/h') * (2, 'day')
        >>> s.has_dimension('m')
        True
        >>> print(s.ounit('cm'))
        576000000.0000 cm
        """
        o = self.sunit2mag(u)
        return self.unit == o.unit

    def _mult_by(self, m: Any) -> None:
        m = self.coerce(m)
        self.val *= m.val
        for i in range(len(self.unit)):
            self.unit[i] = self.unit[i] + m.unit[i]
        self.out_unit = None

    def _div_by(self, m: Any) -> None:
        m = self.coerce(m)
        self.val /= m.val
        for i in range(len(self.unit)):
            self.unit[i] = self.unit[i] - m.unit[i]
        self.out_unit = None

    def ounit(self, unit: Optional[str] = None) -> "Magnitude":
        """Set the preferred unit for output, returning the Magnitude.

        >>> a = mg(1, 'kg m2 / s2')
        >>> print(a)
        1.0000 kg m2 / s2
        >>> print(a.ounit('J'))
        1.0000 J
        >>> print(a)
        1.0000 J
        """
        self.out_unit = unit
        if unit:
            self.out_factor = self.sunit2mag(unit)
            if self.out_factor.unit != self.unit:
                raise IncompatibleUnitsError(
                    self.out_factor.unit, self.unit, "setting output unit"
                )
        else:
            self.out_factor = None
        return self

    def to_base_units(self) -> "Magnitude":
        """Forgets about the output unit and goes back to base units:

        >>> a = mg(10, 'km')
        >>> print(a)
        10.0000 km
        >>> print(a.to_base_units())
        10000.0000 m
        """
        self.out_unit = None
        self.out_factor = None
        return self

    def output_prec(self, prec: int) -> "Magnitude":
        """Set the output precision for the Magnitude.

        If not set, the the module's default will be used, set and
        queried with output_precision(prec).

        >>> a = mg(5, 'm3') ** (1/3.)  # Careful with precedence of **
        >>> print(a)
        1.7100 m
        >>> print(a.output_prec(1))
        1.7 m
        """
        self.oprec = prec
        return self

    def output_format(self, oformat: str) -> "Magnitude":
        """Set the output format for the Magnitude.

        If not set, the module's default will be used, set and queried
        with default_format(fmt).  Default value is "%.*f".  The star
        will be replaced by the expected output precision.

        >>> a = mg(5, 'm2').sqrt()
        >>> print(a)
        2.2361 m
        >>> print(a.output_format("%03d"))
        002 m
        """
        self.oformat = oformat
        return self

    def coerce(self, m: Any) -> "Magnitude":
        """Force tuples or numbers into Magnitude."""
        if isinstance(m, Magnitude):
            return m

        if isinstance(m, tuple):
            if len(m) == 2:
                r = Magnitude(m[0])
                r._mult_by(self.sunit2mag(m[1]))
                return r
            elif len(m) == 1:
                return Magnitude(m[0])
            else:
                raise ConversionError(f"Cannot coerce tuple of length {len(m)} to Magnitude")
        elif isinstance(m, numbers.Number):
            return Magnitude(m)
        else:
            raise ConversionError(f"Cannot coerce {type(m).__name__} to Magnitude")

    def __add__(self, m: Any) -> "Magnitude":
        """Add Magnitude instances.

        >>> print(mg(10, 'm') + (20, 'km') + (30, 'lightyear'))
        283821914177444000.0000 m
        """
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "addition")
        r = self.copy()
        r.val += m.val
        return r

    def __radd__(self, m):
        """Add Magnitude instances.  See __add__."""
        return self.__add__(m)

    def __iadd__(self, m):
        """Add Magnitude instances.  See __add__."""
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "addition")
        self.val += m.val
        return self

    def __sub__(self, m):
        """Substract Magnitude instances.

        >>> print(mg(20, 'm/s') - (1, 'km/h'))
        19.7222 m / s
        """
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "subtraction")
        r = self.copy()
        r.val -= m.val
        return r

    def __rsub__(self, m):
        """Substract Magnitude instances.  See __sub__."""
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "subtraction")
        r = m.copy()
        r.val -= self.val
        return r

    def __isub__(self, m):
        """Substract Magnitude instances.  See __sub__."""
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "subtraction")
        self.val -= m.val
        return self

    def __mul__(self, m):
        """Multiply Magnitude instances.

        >>> print(mg(10, 'm/s') * (10, 's'))
        100.0000 m
        """
        r = self.copy()
        r._mult_by(m)
        return r

    def __rmul__(self, m):
        """Multiply Magnitude instances.  See __mul__."""
        r = self.copy()
        r._mult_by(m)
        return r

    def __imul__(self, m):
        """Multiply Magnitude instances.  See __mul__."""
        self._mult_by(m)
        return self

    def __div__(self, m):
        """Divide Magnitude instances.

        >>> print(mg(100, 'V') / (10, 'kohm'))
        0.0100 A
        """
        r = self.copy()
        r._div_by(m)
        return r

    def __truediv__(self, m):
        """Divide Magnitude instances when "from __future__ import division"
        is in effect.

        >>> print(mg(100, 'V') / (1, 'kohm'))
        0.1000 A
        """
        r = self.copy()
        r._div_by(m)
        return r

    def __rdiv__(self, m):
        """Divide Magnitude instances.  See __div__."""
        m = self.coerce(m)
        m._div_by(self)
        return m

    def __rtruediv__(self, m):
        """Divide Magnitude instances.  See __div__."""
        m = self.coerce(m)
        m._div_by(self)
        return m

    def __idiv__(self, m):
        """Divide Magnitude instances.  See __div__."""
        self._div_by(m)
        return self

    def __itruediv__(self, m):
        """Divide Magnitude instances.  See __div__."""
        self._div_by(m)
        return self

    def __mod__(self, n):
        """Modulus of a Magnitude by a number or a Magnitude.

        Unit is that of the left hand side operator.

        >>> print(mg(10, 'm/s') % 3)
        1.0000 m / s
        >>> print(mg(10, 'm/s') % (3, 'W'))
        1.0000 m / s
        """
        r = self.copy()
        r.val = r.val % ensmg(n).toval()
        return r

    def __imod__(self, n):
        """Modulus of a Magnitude by a number or a Magnitude.  See __mod__."""
        self.val %= n.val
        for i in range(len(self.unit)):
            self.unit[i] = self.unit[i] - n.unit[i]
        self.out_unit = None
        return self

    def __rmod__(self, m):
        """Reverse modulus operation."""
        m = self.coerce(m)
        r = m.copy()
        r.val = r.val % self.toval()
        return r

    def __floordiv__(self, m):
        """Floordiv of two Magnitude instances.

        >>> print(mg(10, 'm/s') // (3, 's'))
        3.0000 m / s2
        >>> print(mg(-10, 'm/s') // (3, 'm'))
        -4.0000 1 / s
        """
        r = self.copy()
        r._div_by(m)
        r.val = math.floor(r.val)
        return r

    def __ifloordiv__(self, m):
        """Floordiv of two Magnitude instances. See __floordiv__."""
        self._div_by(m)
        self.val = math.floor(self.val)
        return self

    def __rfloordiv__(self, m):
        """Reverse floordiv of two Magnitude instances."""
        m = self.coerce(m)
        m._div_by(self)
        m.val = math.floor(m.val)
        return m

    def __divmod__(self, m):
        """Floordiv and remainder of two Magnitude instances.

        >>> [ str(i) for i in divmod(mg(10, 'm/s'), (3, 's')) ]
        ['3.0000 m / s2', '1.0000 m / s']
        """
        return (self.__floordiv__(m), self.__mod__(m))

    def __rdivmod__(self, m):
        """Floordiv and remainder of two Magnitude instances.
        See __divmod___.
        """
        m = self.coerce(m)
        return (m.__floordiv__(self), m.__mod__(self))

    def __pow__(self, n, modulo=None):
        """Return a Magnitude to the power n.

        If modulo is present return the result modulo it.

        >>> print(mg(10, 'm/s') ** 2)
        100.0000 m2 / s2
        >>> print(pow(mg(10, 'km/h'), mg(2))) # Exponent cannot have dimension
        7.7160 m2 / s2
        >>> print(pow(mg(10, 'm/s'), 2, 3))
        1.0000 m2 / s2
        """
        r = self.copy()
        if modulo and (r.val == math.floor(r.val)):  # it's an integer
            # might have been converted to float during creation,
            # modulo only works when all are int
            r.val = int(r.val)
        if isinstance(n, Magnitude):  # happens when called as a ** n
            if not n.dimensionless():
                raise MagnitudeError(f"Cannot use a dimensional number as exponent, {n}")
            n = n.val
        r.val = pow(r.val, n, modulo)
        for i in range(len(r.unit)):
            r.unit[i] *= n
        return r

    def __ipow__(self, n, modulo=None):
        """Power of a Magnitude.  See __pow___."""
        if modulo and (self.val == math.floor(self.val)):  # it's an integer
            # might have been converted to float during creation,
            # modulo only works when all are int
            self.val = int(self.val)
        if isinstance(n, Magnitude):
            if not n.dimensionless():
                raise MagnitudeError(f"Cannot use a dimensional number as exponent, {n}")
            n = n.val
        self.val = pow(self.val, n, modulo)
        for i in range(len(self.unit)):
            self.unit[i] *= n
        return self

    def __neg__(self):
        """Multiply by -1 the value of the Magnitude."""
        r = self.copy()
        r.val = -r.val
        return r

    def __pos__(self):
        """Unary plus operator."""
        return self.copy()

    def __abs__(self):
        """Absolute value of a Magnitude.

        >>> print(abs(mg(-10, 'm')))
        10.0000 m
        """
        r = self.copy()
        r.val = abs(r.val)
        return r

    def __eq__(self, m):
        """Check whether two Magnitude instances with the same dimensions are equal.

        >>> print(mg(1, 'km') == (1000, 'm'))
        True
        """
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "comparison")
        return self.val == m.val

    def __ne__(self, m):
        """Check whether two Magnitude instances with the same dimensions are
        not equal.

        >>> print(mg(1, 'km') != (1000, 'm'))
        False
        """
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "comparison")
        return self.val != m.val

    def __gt__(self, m):
        """Compare two Magnitude instances with the same dimensions.

        >>> print(mg(10, 'm/s') > (10, 'km/h'))
        True
        """
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "comparison")
        return self.val > m.val

    def __ge__(self, m):
        """Compare two Magnitude instances with the same dimensions.

        >>> print(mg(10, 'm/s') >= (10, 'km/h'))
        True
        """
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "comparison")
        return self.val >= m.val

    def __lt__(self, m):
        """Compare two Magnitude instances with the same dimensions.

        >>> print(mg(9, 'm/s') < (10, 'km/h'))
        False
        """
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "comparison")
        return self.val < m.val

    def __le__(self, m):
        """Compare two Magnitude instances with the same dimensions.

        >>> print(mg(10, 'm/s') <= (10, 'km/h'))
        False
        """
        m = self.coerce(m)
        if m.unit != self.unit:
            raise IncompatibleUnitsError(m.unit, self.unit, "comparison")
        return self.val <= m.val

    def __int__(self):
        """Return the value of a Magnitude coerced to integer.

        Note that this will happen to the value in the default output unit:

        >>> print(int(mg(10.5, 'm/s')))
        10
        >>> print(int(mg(10.5, 'm/s').ounit('km/h')))
        37
        """
        return int(self.toval())

    def __long__(self):
        """Return the value of a Magnitude coerced to long.  See __int__."""
        return int(self.toval())

    def __float__(self):
        """Return the value of a Magnitude coerced to float.  See __int__."""
        return float(self.toval())

    def ceiling(self):
        """Ceiling of a Magnitude's value in canonical units.

        >>> print(mg(10.2, 'm/s').ceiling())
        11.0000 m / s
        >>> print(mg(3.6, 'm/s').ounit('km/h').ceiling())
        4.0000 m / s
        >>> print(mg(50.3, 'km/h').ceiling())
        14.0000 m / s
        """
        r = self.copy(with_format=False)
        r.val = math.ceil(r.val)
        return r

    def floor(self):
        """Floor of a Magnitude's value in canonical units.

        >>> print(mg(10.2, 'm/s').floor())
        10.0000 m / s
        >>> print(mg(3.6, 'm/s').ounit('km/h').floor())
        3.0000 m / s
        >>> print(mg(50.3, 'km/h').floor())
        13.0000 m / s
        """
        r = self.copy()
        r.val = math.floor(r.val)
        return r

    def round(self):
        """Round a Magnitude's value in canonical units.

        >>> print(mg(10.2, 'm/s').round())
        10.0000 m / s
        >>> print(mg(3.6, 'm/s').ounit('km/h').round())
        4.0000 m / s
        >>> print(mg(50.3, 'km/h').round())
        14.0000 m / s
        """
        r = self.copy()
        r.val = round(r.val)
        return r

    def to_bits(self):
        return Magnitude(math.ceil(math.log(self.val) / math.log(2.0)), b=1)

    def sqrt(self):
        """Square root of a magnitude.

        >>> print(mg(4, 'm2/s2').sqrt())
        2.0000 m / s
        >>> print(mg(2, 'm/s').sqrt())
        1.4142 m0.5 / s0.5
        """
        return self**0.5


# Some helper functions


def mg(v: Any, unit: str = "", ounit: str = "") -> Magnitude:
    """Builds a Magnitude from a number and a units string.  Specify
    the preferred output unit with ounit (by default equals to unit).

    >>> print(mg(10, 'm/s'))
    10.0000 m/s
    >>> a = mg(10, 'm/s', 'km/h')
    >>> print(a)
    36.0000 km/h
    >>> a = mg(1, 'B')
    >>> print(a)
    1.0000 B
    >>> print(a.ounit('b'))
    8.0000 b
    >>> a = mg(1024, 'B')
    >>> print(a.ounit('b'))
    8192.0000 b
    >>> print(a.ounit('KiB'))
    1.0000 KiB
    """
    m = Magnitude(v)
    if unit:
        u = m.sunit2mag(unit)
        m._mult_by(u)
    if not ounit:
        ounit = unit
    return m.ounit(ounit)


def ensmg(
    m: Optional[Union[float, tuple[float], tuple[float, str], Magnitude]],
    unit: str = "",
) -> Magnitude:
    """Converts something to a Magnitude.

    >>> print(ensmg(10, 'Hz'))
    10.0000 Hz
    >>> print(ensmg(ensmg(1000, 'Hz')))
    1000.0000 Hz
    >>> a = (4, 'mol')
    >>> print(ensmg(a))
    4.0000 mol
    >>> a = mg(1024, 'Pa')
    >>> print(ensmg(a))
    1024.0000 Pa
    >>> f = ensmg((10, 'Pa')) * (10, 'm2')
    >>> print(f.ounit('N'))
    100.0000 N
    """
    if not isinstance(m, Magnitude):
        if isinstance(m, tuple):
            if len(m) == 2:
                return mg(m[0], m[1], unit)
            elif (len(m) == 1) and isinstance(m[0], numbers.Number):
                if unit:
                    return mg(m[0], unit)
                return Magnitude(m[0])
            else:
                raise ConversionError(f"Can't convert {m} to Magnitude")
        elif isinstance(m, numbers.Number):
            if unit:
                return mg(m, unit)
            return Magnitude(m)
        else:
            raise ConversionError(f"Can't convert {m} to Magnitude")
    else:
        return m


def new_mag(indicator: str, mag: Magnitude) -> None:
    """Define a new magnitude understood by the package.

    Defines a new magnitude type by giving it a name (indicator) and
    its equivalence in the form of an already understood magnitude.

    >>> new_mag('mile', mg(160934.4, 'cm'))
    >>> print(mg(100, 'mile/h').ounit('km/h'))
    160.9344 km/h
    """
    _mags[indicator] = mag


# Finally, define the Magnitudes and initialize _mags.


def _init_mags():
    # Magnitudes for the base SI units
    new_mag("m", Magnitude(1.0, m=1))
    new_mag("s", Magnitude(1.0, s=1))
    new_mag("K", Magnitude(1.0, K=1))
    new_mag("kg", Magnitude(1.0, kg=1))
    new_mag("A", Magnitude(1.0, A=1))
    new_mag("mol", Magnitude(1.0, mol=1))
    new_mag("cd", Magnitude(1.0, cd=1))
    new_mag("$", Magnitude(1.0, dollar=1))
    new_mag("dollar", Magnitude(1.0, dollar=1))
    new_mag("b", Magnitude(1.0, b=1))  # bit

    # Magnitudes for derived SI units
    new_mag("B", Magnitude(8.0, b=1))
    new_mag("rad", Magnitude(1.0))  # radian
    new_mag("sr", Magnitude(1.0))  # steradian
    new_mag("Hz", Magnitude(1.0, s=-1))  # hertz
    new_mag("g", Magnitude(1e-3, kg=1))  # gram
    new_mag("N", Magnitude(1.0, m=1, kg=1, s=-2))  # newton
    new_mag("Pa", Magnitude(1.0, m=-1, kg=1, s=-2))  # pascal
    new_mag("J", Magnitude(1.0, m=2, kg=1, s=-2))  # joule
    new_mag("W", Magnitude(1.0, m=2, kg=1, s=-3))  # watt
    new_mag("C", Magnitude(1.0, s=1, A=1))  # coulomb
    new_mag("V", Magnitude(1.0, m=2, kg=1, s=-3, A=-1))  # volt
    new_mag("F", Magnitude(1.0, m=-2, kg=-1, s=4, A=2))  # farad, C/V
    new_mag("ohm", Magnitude(1.0, m=2, kg=1, s=-3, A=-2))  # ohm, V/A
    new_mag("S", Magnitude(1.0, m=-2, kg=-1, s=3, A=2))  # siemens, A/V, el cond
    new_mag("Wb", Magnitude(1.0, m=2, kg=1, s=-2, A=-1))  # weber, V.s, mag flux
    new_mag("T", Magnitude(1.0, kg=1, s=-2, A=-1))  # tesla, Wb/m2, mg flux dens
    new_mag("H", Magnitude(1.0, m=2, kg=1, s=-2, A=-2))  # henry, Wb/A, induct.
    new_mag("degC", Magnitude(1.0, K=1))  # celsius, !!
    new_mag("lm", Magnitude(1.0, cd=1))  # lumen, cd.sr (=cd)), luminous flux
    new_mag("lux", Magnitude(1.0, m=-2, cd=1))  # lux, lm/m2, illuminance
    new_mag("Bq", Magnitude(1.0, s=-1))  # becquerel, activity of a radionulide
    new_mag("Gy", Magnitude(1.0, m=2, s=-2))  # gray, J/kg, absorbed dose
    new_mag("Sv", Magnitude(1.0, m=2, s=-2))  # sievert, J/kg, dose equivalent
    new_mag("kat", Magnitude(1.0, s=-1, mol=1))  # katal, catalitic activity

    # length
    new_mag("'", Magnitude(0.3048, m=1))  # feet
    new_mag("ft", Magnitude(0.3048, m=1))  # feet
    new_mag("inch", Magnitude(0.0254, m=1))  # inch
    new_mag('"', Magnitude(0.0254, m=1))  # inch
    new_mag("lightyear", Magnitude(2.99792458e8 * 365.25 * 86400, m=1))

    # volume
    new_mag("l", Magnitude(0.001, m=3))

    # time
    # year is tropical year, "the mean interval between vernal
    # equinoxes.  Differs from the sidereal year by 1 part in 26000
    # due to precession of the earth about its rotational axis
    # combined with precession of the perihelion of the earth's orbit"
    # (from units.dat).
    new_mag("year", Magnitude(31556925.974678401, s=1))
    new_mag("day", Magnitude(86400, s=1))
    new_mag("h", Magnitude(3600, s=1))
    new_mag("min", Magnitude(60, s=1))

    # Resolution
    new_mag("dpi", Magnitude(1.0 / 0.0254, m=-1))
    new_mag("lpi", Magnitude(1.0 / 0.0254, m=-1))

    # Velocity
    new_mag("ips", Magnitude(0.0254, m=1, s=-1))
    new_mag("c", Magnitude(2.99792458e8, m=1, s=-1))

    # Acceleration
    new_mag("gravity", Magnitude(9.80665, m=1, s=-2))

    # Coverage
    new_mag("gsm", Magnitude(0.001, kg=1, m=-2))


if not _mags:
    _init_mags()


if __name__ == "__main__":
    import doctest

    doctest.testmod()
