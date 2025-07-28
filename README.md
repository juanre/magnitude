# Magnitude

A Python library for computing with numbers with units.

## About

A magnitude is a number with a unit, like 10 km/h. Units can be any of the SI units, plus a bunch of non-SI, bits, dollars, and any combination of them. They can include the standard SI prefixes. Magnitude can operate with physical quantities, parse their units, and print them. You don't have to worry about unit consistency or conversions; everything is handled transparently. By default output is done in basic SI units, but you can specify any output unit, as long as it can be reduced to the basic units of the physical quantity.

## Installation

```bash
uv add magnitude # recommended

# or with pip

pip install magnitude
```

## Quick Example

```python
from magnitude import mg

# Basic usage
print("10 m/s ** 2 ->", mg(10, 'm/s') ** 2)
# 10 m/s ** 2 -> 100.0000 m2 / s2

# Unit conversions
speed = mg(10, 'm/s', 'km/h')
speed
# 36.0000 km/h

# Dimensional analysis
time_squared = mg(10, 'm') * 2 / (10, 'm/s2')
time_squared
# 2.0000 s2
back_to_time = time_squared ** 0.5
back_to_time
# 1.4142 s

# Natural constants
year = mg(1, "lightyear") / (1, "c")
year.ounit("year")
# 1.0000 year
year.ounit('day')
# 365.2500 day
```

## Features

- Full support for SI units and prefixes
- Automatic unit conversion and dimensional analysis
- Support for derived units (joules, watts, etc.)
- Binary prefixes (Ki, Mi, Gi, etc.)
- Currency units
- Unicode superscript support for unit notation (m², s⁻¹, etc.)
- Extensible - add your own units
- Pure Python, no dependencies

## Basic Units

The basic units understood by magnitude are:

- `$` - dollar
- `A` - ampere
- `b` - bit
- `cd` - candela
- `K` - degrees Kelvin
- `kg` - kilograms
- `m` - meters
- `mol` - amount of substance
- `s` - seconds

## Derived Units

From the basic units, magnitude supports many derived units:

- `Bq` - becquerel
- `C` - coulomb
- `c` - speed of light (m/s)
- `day` - day
- `degC` - degree Celsius
- `dpi` - dots per inch
- `F` - farad
- `ft` - feet ("'" is also acceptable)
- `g` - gram
- `gravity` - acceleration due to gravity (m/s²)
- `Gy` - gray
- `H` - henry
- `h` - hour
- `Hz` - Hertz
- `inch` - inch ('"' is also acceptable)
- `ips` - inches per second
- `J` - joule
- `kat` - katal
- `l` - liter
- `lightyear` - light year
- `lm` - lumen
- `lpi` - lines per inch
- `lux` - lux
- `min` - minute
- `N` - newton
- `ohm` - ohm
- `Pa` - pascal
- `S` - siemens
- `Sv` - sievert
- `T` - tesla
- `V` - volt
- `W` - watt
- `Wb` - weber
- `year` - year
- `B` - byte

Two magnitudes have no units: `rad` (radian - unit of plane angle) and `sr` (steradian - unit of solid angle).

## Scale Prefixes

Any unit can be augmented with these scale prefixes:

### Decimal Prefixes
- `y` - yocto (10⁻²⁴)
- `z` - zepto (10⁻²¹)
- `a` - atto (10⁻¹⁸)
- `f` - femto (10⁻¹⁵)
- `p` - pico (10⁻¹²)
- `n` - nano (10⁻⁹)
- `u` - micro (10⁻⁶)
- `m` - milli (10⁻³)
- `c` - centi (10⁻²)
- `d` - deci (10⁻¹)
- `k` - kilo (10³)
- `M` - mega (10⁶)
- `G` - giga (10⁹)
- `T` - tera (10¹²)
- `P` - peta (10¹⁵)
- `E` - exa (10¹⁸)
- `Z` - zetta (10²¹)
- `Y` - yotta (10²⁴)

### Binary Prefixes
- `Ki` - Kibi (2¹⁰)
- `Mi` - Mebi (2²⁰)
- `Gi` - Gibi (2³⁰)
- `Ti` - Tebi (2⁴⁰)
- `Pi` - Pebi (2⁵⁰)
- `Ei` - Exbi (2⁶⁰)

## Unicode Superscripts

Magnitude supports Unicode superscripts for both input and output of units:

### Input
You can use Unicode superscripts when creating magnitudes:

```python
from magnitude import mg

# These are equivalent
area1 = mg(10, 'm²')
area2 = mg(10, 'm2')

# Works with negative exponents too
frequency = mg(440, 's⁻¹')  # Same as mg(440, 's-1') or mg(440, 'Hz')

# And multi-digit exponents
special = mg(1, 'kg¹²')  # Same as mg(1, 'kg12')
```

### Output
By default, output uses regular ASCII notation. You can enable Unicode superscripts:

```python
from magnitude import mg, unicode_superscript

m = mg(10, 'm2/s2')
print(m)  # 10.0000 m2/s2

# Enable Unicode superscripts
unicode_superscript(True)
print(m)  # 10.0000 m²/s²

# Works with all units
print(mg(1, 'kg') / mg(1, 'm3'))  # 1.0000 kg / m³
```

### Environment Variable
You can set the default behavior using the `MAGNITUDE_UNICODE_SUPERSCRIPTS` environment variable:

```bash
# Enable Unicode superscripts by default
export MAGNITUDE_UNICODE_SUPERSCRIPTS=1  # or 'true', 'yes', 'on'

# In Python
from magnitude import mg
print(mg(10, 'm2'))  # Will print: 10.0000 m²
```

## Defining New Magnitudes

You can define new magnitudes by instantiating the Magnitude class. For example, to define pounds:

```python
from magnitude import Magnitude, mg, new_mag

# A pound is 0.45359237 kilograms
lb = Magnitude(0.45359237, kg=1)

# Register it in the system
new_mag('lb', lb)

# Now you can use it
me = mg(180, 'lb')
print(me.ounit('kg'))  # 81.6466 kg
```

## API Reference

### Main Classes and Functions

- `Magnitude` - The main class for numbers with units
- `mg(value, unit, ounit='')` - Construct a Magnitude
- `ensmg(m, unit='')` - Convert something to a Magnitude
- `new_mag(indicator, mag)` - Register a new magnitude unit
- `MagnitudeError` - Exception for magnitude errors

### Output Formatting

- `output_precision(prec)` - Set/get output precision (default: 4)
- `output_units(enabled)` - Enable/disable unit output
- `default_format(fmt)` - Set/get default output format

## Documentation

Full documentation is available at: https://juanreyero.com/open/magnitude/

## References

- [NIST Units](http://physics.nist.gov/cuu/Units/units.html)
- [SI Units on Wikipedia](http://en.wikipedia.org/wiki/SI)
- [GNU Units](http://www.gnu.org/software/units/units.html)

This code was inspired by [Novak's units code](http://www.cs.utexas.edu/users/novak/units.html) and its [associated paper](http://www.cs.utexas.edu/users/novak/units95.html).

## License

Copyright (c) 2006-2024 Juan Reyero (https://juanreyero.com).

Licensed under the MIT License. See LICENSE for details.
