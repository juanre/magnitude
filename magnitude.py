# magnitude  -- a module for computing with numbers with units.
#
# Version 0.9.4, November 2009
#
# Copyright (C) 2006-2009 Juan Reyero (http://juanreyero.com).
# 
# Licensed under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied.  See the License for the specific
# language governing permissions and limitations under the
# License.
# 
# Home page: http://juanreyero.com/magnitude/

"""
A physical quantity is a number with a unit, like 10 km/h. Units can be any
of the SI units, plus a bunch of non-SI, bits, dollars, and any combination
of them. They can include the standard SI prefixes. Magnitude can operate
with physical quantities, parse their units, and print them. You don't have
to worry about unit consistency or conversions; everything is handled
transparently. By default output is done in basic SI units, but you can
specify any output unit, as long as it can be reduced to the basic units of
the physical quantity.

The basic units understood by the magnitude module are:

    indicator    meaning
    ---------    -------
    $            dollar ('dollar' is also acceptable)
    A            amperes
    b            binary digit
    cd           candela
    K            degrees Kelvin
    kg           kilograms
    m            meters
    mol          amount of substance
    s            seconds

From these basic units you can derive many other units.  The magnitude
package predefines these derived units:

    Bq           becquerel
    C            coulomb
    c            speed of light (m/s)
    day
    degC         degree Celsius
    dpi          dots per inch
    F            farad
    ft           feet ("'" is also acceptable)
    g            gram
    gravity      acceleration due to gravity (m/s**2)
    Gy           gray
    H            henry
    h            hour
    Hz           Hertz
    inch         ('"' is also acceptable)
    ips          inches per second
    J            joule
    kat          byte
    l            liter
    lightyear    light year
    lm           lumen
    lpi          lines per inch
    lux
    min          minute
    N            newton
    ohm
    Pa           pascal
    S            siemens
    Sv           sievert
    T            tesla
    V            volt
    W            watt
    Wb           weber
    year
    B            byte

Two magnitudes have no units, 'rad' (radian - unit of plane angle) and 'sr'
(steradian - unit of solid angle).

Any of the above units can be augmented with the following set of scale
prefixes:

    letter     scale    name
    ------     -----    ----
    y          1e-24    yocto
    z          1e-21    zepto
    a          1e-18    atto
    f          1e-15    femto
    p          1e-12    pico
    n          1e-9     nano
    u          1e-6     micro
    m          1e-3     mili
    c          1e-2     centi
    d          1e-1     deci
    k          1e3      kilo
    Ki         2^10     Kibi
    M          1e6      mega
    Mi         2^20     Mebi
    G          1e9      giga
    Gi         2^30     Gibi
    T          1e12     tera
    Ti         2^40     Tebi
    P          1e15     peta
    Pi         2^50     Pebi
    E          1e18     exa
    Ei         2^60     Exbi
    Z          1e21     zetta
    Y          1e24     yotta

You can define new magnitudes by instantiating the Magnitude class.  Suppose
you want to define pounds as a magnitude and associate with it the unit
'lb'.  A pound is 0.45359237 kilograms, so we have

    >>> lb = Magnitude(0.45359237, kg=1)

To make it recognized automatically you also have to introduce it to
the system:

    >>> new_mag('lb', lb)

You can then use it as you would any other predefined physical quantity:

    >>> me = mg(180, 'lb')
    >>> print me.ounit('kg').toval()
    81.6466266

The following online references provide more detail about physical units and
the SI system.

    http://physics.nist.gov/cuu/Units/units.html
    http://en.wikipedia.org/wiki/SI
    http://www.gnu.org/software/units/units.html for units.dat
    http://www.cip.physik.uni-muenchen.de/~tf/misc/etools.lisp

This code was very much inspired by
    http://www.cs.utexas.edu/users/novak/units.html
and its associated paper,
    http://www.cs.utexas.edu/users/novak/units95.html

Bits and bytes (2009-11-03)
---------------------------

A previous version of the library used "bit" for bit and "b" for byte,
leaving B for Bel.  Following Michael Scheper's suggestion we follow
now IEEE 1541 and use "b" for bit and "B" for byte.  If the need
arises I'll implement ad-hoc esupport for dB, but for the time being
there is none.
"""

import re, math
import types

# Base magnitude names and prefixes.  The _mags dictionary, initialized
# at the end, will contain all the known magnitudes.  Units are
# 9-element arrays, each element the exponent of the unit named by the
# Uname in the same position. 

class MagnitudeError(Exception):
    pass

_mags = {}
_unames = ['m', 's', 'K', 'kg', 'A', 'mol', 'cd', '$', 'b']
_prefix = {'y': 1e-24,  # yocto
           'z': 1e-21,  # zepto
           'a': 1e-18,  # atto
           'f': 1e-15,  # femto
           'p': 1e-12,  # pico
           'n': 1e-9,   # nano
           'u': 1e-6,   # micro
           'm': 1e-3,   # mili
           'c': 1e-2,   # centi
           'd': 1e-1,   # deci
           'k': 1e3,    # kilo
           'M': 1e6,    # mega
           'G': 1e9,    # giga
           'T': 1e12,   # tera
           'P': 1e15,   # peta
           'E': 1e18,   # exa
           'Z': 1e21,   # zetta
           'Y': 1e24,   # yotta

           # Binary prefixes, approved by the International
           # Electrotechnical Comission in 1998.  Since then, kb means
           # 1000 bytes; for 1024 bytes use Kib (note the capital K in
           # the binary version, and the lower case for the b of byte,
           # see comment in byte definition below).
           'Ki': 2 ** 10, # Kibi (<- kilo, 10^3)
           'Mi': 2 ** 20, # Mebi (<- mega, 10^6)
           'Gi': 2 ** 30, # Gibi (<- giga, 10^9)
           'Ti': 2 ** 40, # Tebi (<- tera, 10^12)
           'Pi': 2 ** 50, # Pebi (<- peta, 10^15)
           'Ei': 2 ** 60  # Exbi (<- exa, 10^18)
           }


###### Default print formatting options

_prn_format = "%.*f"
_prn_prec = 4
_prn_units = True

def default_format(fmt=None):
    """Get or set the default ouptut format.  Include a fmt if and
    where you need to specify the output precission.  Defaults to
    %.*f, where the * stands for the precission.  Returns the default
    format.

    >>> print mg(2, 'm2').sqrt()
    1.4142 m
    >>> default_format("%.2f")
    '%.2f'
    >>> print mg(2, 'm2').sqrt()
    1.41 m
    """
    global _prn_format
    if fmt is not None:
        _prn_format = fmt
    return _prn_format

def output_precision(prec=None):
    """Get or set the output precission.  Default is 4. Returns it.

    >>> default_format("%.*f") # This is the out-of-the-box default
    '%.*f'
    >>> print mg(2, 'm2').sqrt()
    1.4142 m
    >>> output_precision(6)
    6
    >>> print mg(2, 'm2').sqrt()
    1.414214 m
    >>> output_precision(4)
    4
    """
    global _prn_prec
    if prec is not None:
        _prn_prec = prec
    return _prn_prec

def output_units(un=None):
    """Get or set the units output.  Default is True.  If False output
    is only numbers.

    >>> print mg(2, 'day')
    2.0000 day
    >>> output_units(False)
    False
    >>> print mg(2, 'day').ounit('s')
    172800.0000
    """
    global _prn_units
    if un is not None:
        _prn_units = un
    return _prn_units


###### Resolution areas

def res2num(res):
    match = re.search(r'(\d+)x(\d+)', res)
    if match:
        return int(match.group(1)), int(match.group(2))
    if (res[0] == '[') and (res[-1] == ']'):
        return (int(res[1:-1]), int(res[1:-1]))

def _isres(res):
    return (len(res) > 2) and (res[0] == '[') and (res[-1] == ']')

def res2m2(res):
    """Bracketed resolutions are used in the printing industry, to
    denote the area of a pixel.  Can be like [300x1200] or like [600]
    (=[600x600]), meaning the area of square pixels of size 1"/300 x
    1"/1200 and 1"/600 x 1"/600.  The square brackes are intended to
    show that we are talking about areas.  This function converts them
    to square meters.

    >>> res2m2("[600x600]")
    1.7921111111111111e-09
    >>> res2m2("[600]")
    1.7921111111111111e-09
    >>> res2m2("[150x300]")
    1.4336888888888889e-08
    """
    hr, vr = res2num(res)
    return 0.0254 * 0.0254 / (vr * hr)


# Definition of the magnitude type.  Includes operator overloads.

def numberp(n):  ## Python has to have a decent way to do this!
    return (isinstance(n, types.ComplexType) or
            isinstance(n, types.FloatType) or
            isinstance(n, types.IntType) or
            isinstance(n, types.LongType))

class Magnitude:
    def __init__(self, val, m=0, s=0, K=0, kg=0, A=0, mol=0, cd=0, dollar=0,
                 b=0):
        self.val = val
        self.unit = [m, s, K, kg, A, mol, cd, dollar, b]
        self.out_unit = None
        self.out_factor = None
        self.oprec = None
        self.oformat = None

    def copy(self):
        """Builds and returns a copy of a magnitude, except the
        default output unit, output factor, output precision and
        output format.
        >>> a = mg(1000, 'um')
        >>> print a
        1000.0000 um
        >>> b = a.copy()
        >>> print b
        0.0010 m
        """
        return Magnitude(self.val, *self.unit)

    def toval(self, ounit=''):
        """Returns the numeric value of a magnitude in ounit or in its
        default output unit.
        >>> v = mg(100, 'km/h')
        >>> v.toval()
        100.0
        >>> v.toval(ounit='m/s')
        27.777777777777779
        """
        m = self.copy()
        if not ounit:
            ounit = self.out_unit
        if ounit:
            out_factor = self.sunit2mag(ounit)
            m.div_by(out_factor)
        return m.val

    def __str__(self):
        oformat = self.oformat
        oprec = self.oprec
        if oprec is None:
            oprec = _prn_prec
        if oformat is None:
            oformat = _prn_format
        if self.out_unit:
            m = self.copy()
            m.div_by(self.out_factor)
            if '*' in oformat:  # requires the precision arg
                st = oformat % (oprec, m.val)
            else:
                st = oformat % (m.val)
            if _prn_units:
                return st + ' ' + self.out_unit.strip()
            return st

        if '*' in oformat:
            st = oformat % (oprec, self.val)
        else:
            st = oformat % (self.val)
            
        if not _prn_units:
            return st

        u = self.unit
        num = ' '  # numerator
        for i in range(len(_unames)):
            if u[i] > 1:
                num = num + _unames[i] + str(u[i]) + ' '
            elif u[i] == 1:
                num = num + _unames[i] + ' '
        den = ''  # denominator
        for i in range(len(_unames)):
            if u[i] < -1:
                den = den + _unames[i] + str(-u[i]) + ' '
            elif u[i] == -1:
                den = den + _unames[i] + ' '
        if den:
            if num == ' ':
                num += '1 '
            st += (num + '/ ' + den)
        elif num != ' ':
            st += num        
        return st.strip()

    def term2mag(self, s):
        """Converts a string with units to a Magnitude.  Can't divide: use
        with the numerator and the denominator separately (hence the term).

        """
        m = Magnitude(1.0)
        units = re.split(r'\s', s)
        for u in units:
            if re.search(r'[^\s]', u):
                exp = 1
                if re.search(r'\d$', u):
                    exp = int(u[-1])
                    u = u[0:-1]
                if _mags.has_key(u):
                    u = _mags[u].copy()
                elif ((len(u)>=3) and _prefix.has_key(u[0:2]) and
                      _mags.has_key(u[2:])):
                    pr = _prefix[u[0:2]]
                    u = _mags[u[2:]].copy();  u.val = pr * u.val
                elif ((len(u)>=2) and _prefix.has_key(u[0]) and
                      _mags.has_key(u[1:])):
                    pr = _prefix[u[0]]
                    u = _mags[u[1:]].copy();  u.val = pr * u.val
                elif _isres(u):
                    u = Magnitude(res2m2(u), m=2)
                elif u == '':
                    u = Magnitude(1.0)
                else:
                    raise MagnitudeError("Don't know about unit %s" % u)
                for i in range(exp):
                    m.mult_by(u)
        return m

    def sunit2mag(self, unit=''):
        """This is the good one: uses term2mag to convert a string
        with units, possibly including a / to separate a numerator and
        a denominator, to a Magnitude.

        >>> a = mg(1, '')
        >>> a.sunit2mag('m/s').toval()
        1.0
        >>> a.sunit2mag('km/h').toval()
        0.27777777777777779
        >>> print a.sunit2mag('W h')
        3600.0000 m2 kg / s2
        >>> print a.sunit2mag('W h').ounit('J')
        3600.0000 J
        """
        m = Magnitude(1.0)
        if unit:
            q = re.split(r'/', unit)
            if re.search(r'[^\s]', q[0]):
                m.mult_by(self.term2mag(q[0]))
            if (len(q) == 2) and re.search(r'[^\s]', q[1]):
                m.div_by(self.term2mag(q[1]))
        return m

    def dimensionless(self):
        """True if the magnitude's dimension exponents are all zero. 

        >>> mg(2, 'K').dimensionless()
        False
        >>> mg(2, 'rad').dimensionless()
        True
        """
        return self.unit == [0] * 9

    def dimension(self):
        """Return the dimension of the unit. 

        >>> mg(2, 'J').dimension()
        [2, -2, 0, 1, 0, 0, 0, 0, 0]
        """
        return self.unit[:]
    
    def has_dimension(self, u):
        """Returns true if the dimension of the magnitude matches u:

        >>> s = mg(120, 'km/h') * (2, 'day')
        >>> s.has_dimension('m')
        True
        >>> print s.ounit('cm')
        576000000.0000 cm
        """
        o = self.sunit2mag(u)
        return (self.unit == o.unit)

    def mult_by(self, m):
        self.val *= m.val
        for i in range(len(self.unit)):
            self.unit[i] = self.unit[i] + m.unit[i]
        self.out_unit = None

    def div_by(self, m):
        self.val /= m.val
        for i in range(len(self.unit)):
            self.unit[i] = self.unit[i] - m.unit[i]
        self.out_unit = None

    def ounit(self, unit):
        self.out_unit = unit
        self.out_factor = self.sunit2mag(unit)
        if self.out_factor.unit != self.unit:
            raise MagnitudeError("Inconsistent Magnitude units: %s, %s" % 
                                (self.out_factor.unit, self.unit))
        return self

    def output_prec(self, prec):
        self.oprec = prec

    def output_format(self, oformat):
        self.oformat = oformat

    def __coerce__(self, m):
        if not isinstance(m, Magnitude):
            if type(m) == tuple:
                if len(m) == 2:
                    r = Magnitude(m[0])
                    r.mult_by(self.sunit2mag(m[1]))
                    return self, r
                elif len(m) == 1:
                    return self, Magnitude(m[0])
                else:
                    return None
            elif numberp(m):
                return self, Magnitude(m)
            else:
                return None
        else:
            return self, m

    def __add__(self, m):
        if m.unit != self.unit:
            raise MagnitudeError("Incompatible units: %s and %s" %
                                 (m.unit, self.unit))
        r = self.copy()
        r.val += m.val
        return r

    def __radd__(self, m):
        return self.__add__(m)

    def __iadd__(self, m):
        if m.unit != self.unit:
            raise MagnitudeError("Incompatible units: %s and %s" %
                                 (m.unit, self.unit))
        self.val += m.val
        return self
    
    def __sub__(self, m):
        if m.unit != self.unit:
            raise MagnitudeError("Incompatible units: %s and %s" %
                                 (m.unit, self.unit))
        r = self.copy()
        r.val -= m.val
        return r

    def __rsub__(self, m):
        return m.__sub__(self)

    def __isub__(self, m):
        if m.unit != self.unit:
            raise MagnitudeError("Incompatible units: %s and %s" %
                                 (m.unit, self.unit))
        self.val -= m.val
        return self

    def __mul__(self, m):
        r = self.copy()
        r.mult_by(m)
        return r

    def __rmul__(self, m):
        r = self.copy()
        r.mult_by(m)
        return r

    def __imul__(self, m):
        self.mult_by(m)
        return self

    def __div__(self, m):
        r = self.copy()
        r.div_by(m)
        return r

    def __rdiv__(self, m):
        r = self.copy()
        m.div_by(r)
        return m

    def __idiv__(self, m):
        self.div_by(m)
        return self

    def __mod__(self, n):
        r = self.copy()
        r.val = r.val % n.val
        for i in range(len(r.unit)):
            r.unit[i] = r.unit[i] - n.unit[i]
        self.out_unit = None
        return r

    def __imod__(self, n):
        self.val %= n.val
        for i in range(len(self.unit)):
            self.unit[i] = self.unit[i] - n.unit[i]
        self.out_unit = None
        return self

    def __floordiv__(self, m):
        r = self.copy()
        r.div_by(m)
        r.val = math.floor(r.val)
        return r

    def __ifloordiv__(self, m):
        self.div_by(m)
        self.val = math.floor(self.val)
        return self
        
    def __divmod__(self, m):
        return (self.__floordiv__(m), self.__mod__(m))

    def __rdivmod__(self, m):
        return (m.__floordiv__(self), m.__mod__(self))

    def __pow__(self, n, modulo=None):
        r = self.copy()
        if modulo and (r.val == math.floor(r.val)):  # it's an integer
            # might have been converted to float during creation,
            # modulo only works when all are int
            r.val = int(r.val)
        if isinstance(n, Magnitude):  # happens when called as a ** n
            if not n.dimensionless():
                raise MagnitudeError("Cannot use a dimensional number as"
                                     "exponent, %s" % (n))
            n = n.val
        r.val = pow(r.val, n, modulo)
        for i in range(len(r.unit)):
            r.unit[i] *= n
        return r

    def __ipow__(self, n):
        if not n.dimensionless():
            raise MagnitudeError("Cannot use a dimensional number as"
                                 "exponent, %s" % (n))
        n = n.val
        self.val = pow(self.val, n)
        for i in range(len(self.unit)):
            self.unit[i] *= n
        return self

    def __neg__(self):
        r = self.copy()
        r.val = -r.val
        return r

    def __pos__(self):
        return self.copy()

    def __abs__(self):
        r = self.copy()
        r.val = abs(r.val)
        return r

    def __cmp__(self, m):
        if m.unit != self.unit:
            raise MagnitudeError("Incompatible units in comparison: %s and %s" %
                                 (m.unit, self.unit))
        return cmp(self.val, m.val)

    def __int__(self):
        return int(self.toval())

    def __long__(self):
        return long(self.toval())

    def __float__(self):
        return float(self.toval())

    def ceiling(self):
        r = self.copy()
        r.val = math.ceil(r.val)
        return r

    def floor(self):
        r = self.copy()
        r.val = math.floor(r.val)
        return r

    def round(self):
        r = self.copy()
        r.val = round(r.val)
        return r

    def to_bits(self):
        return Magnitude(math.ceil(math.log(self.val) / math.log(2.0)),
                         b=1)

    def sqrt(self):
        return self ** 0.5
        
    
# Some helper functions

def mg(v, unit='', ounit=''):
    """Builds a Magnitude from a number and a units string"""
    m = Magnitude(v)
    if unit:
        u = m.sunit2mag(unit)
        m.mult_by(u)
    if not ounit:
        ounit = unit
    m.ounit(ounit)
    return m

def ensmg(m, unit=''):
    """Ensures that something is a Magnitude."""
    if not isinstance(m, Magnitude):
        if type(m) == tuple:
            if len(m) == 2:
                return mg(m[0], m[1], unit)
            elif (len(m) == 1) and numberp(m[0]):
                if unit:
                    return mg(m[0], unit)
                return Magnitude(m[0])
            else:
                raise MagnitudeError("Can't convert %s to Magnitude" %
                                     (m,))
        elif numberp(m):
            if unit:
                return mg(m, unit)
            return Magnitude(m)
        else:
            raise MagnitudeError("Can't convert %s to Magnitude" %
                                 (m,))
    else:
        return m

# These don't really help much, as it's much easier to use the
# overriden * and / operators.

def mul(m1, *rest):
    m = ensmg(m1)
    for m2 in rest:  m.mult_by(ensmg(m2))
    return m

def div(m1, *rest):
    if rest:
        m = ensmg(m1)
        for m2 in rest:  m.div_by(ensmg(m2))
        return m
    else:
        m = Magnitude(1.0)
        m.div_by(ensmg(m1))
        return m

def new_mag(indicator, mag):
    """Define a new magnitude understood by the package."""
    _mags[indicator] = mag

# Finally, define the Magnitudes and initialize _mags.

def _init_mags():
    # Magnitudes for the base SI units
    new_mag('m', Magnitude(1.0, m=1))
    new_mag('s', Magnitude(1.0, s=1))
    new_mag('K', Magnitude(1.0, K=1))
    new_mag('kg', Magnitude(1.0, kg=1))
    new_mag('A', Magnitude(1.0, A=1))
    new_mag('mol', Magnitude(1.0, mol=1))
    new_mag('cd', Magnitude(1.0, cd=1))
    new_mag('$', Magnitude(1.0, dollar=1))
    new_mag('dollar', Magnitude(1.0, dollar=1))
    new_mag('b', Magnitude(1.0, b=1))

    # Magnitudes for derived SI units
    new_mag('B', Magnitude(8.0, b=1))
    new_mag('rad', Magnitude(1.0))  # radian
    new_mag('sr', Magnitude(1.0))  # steradian
    new_mag('Hz', Magnitude(1.0, s=-1))  # hertz
    new_mag('g', Magnitude(1e-3, kg=1))  # gram
    new_mag('N', Magnitude(1.0, m=1, kg=1, s=-2))  # newton
    new_mag('Pa', Magnitude(1.0, m=-1, kg=1, s=-2))  # pascal    
    new_mag('J', Magnitude(1.0, m=2, kg=1, s=-2))  # joule
    new_mag('W', Magnitude(1.0, m=2, kg=1, s=-3))  # watt
    new_mag('C', Magnitude(1.0, s=1, A=1))  # coulomb
    new_mag('V', Magnitude(1.0, m=2, kg=1, s=-3, A=-1))  # volt
    new_mag('F', Magnitude(1.0, m=-2, kg=-1, s=4, A=2))  # farad, C/V
    new_mag('ohm', Magnitude(1.0, m=2, kg=1, s=-3, A=-2))  # ohm, V/A
    new_mag('S', Magnitude(1.0, m=-2, kg=-1, s=3, A=2))  # siemens, A/V, el cond.
    new_mag('Wb', Magnitude(1.0, m=2, kg=1, s=-2, A=-1))  # weber, V.s, mag. flux
    new_mag('T', Magnitude(1.0, kg=1, s=-2, A=-1))  # tesla, Wb/m2, mg flux dens.
    new_mag('H', Magnitude(1.0, m=2, kg=1, s=-2, A=-2))  # henry, Wb/A, induct.
    new_mag('degC', Magnitude(1.0, K=1))  # celsius, !!
    new_mag('lm', Magnitude(1.0, cd=1))  # lumen, cd.sr (=cd)), luminous flux
    new_mag('lux', Magnitude(1.0, m=-2, cd=1))  # lux, lm/m2, illuminance
    new_mag('Bq', Magnitude(1.0, s=-1))  # becquerel, activity of a radionulide
    new_mag('Gy', Magnitude(1.0, m=2, s=-2))  # gray, J/kg, absorbed dose
    new_mag('Sv', Magnitude(1.0, m=2, s=-2))  # sievert, J/kg, dose equivalent
    new_mag('kat', Magnitude(1.0, s=-1, mol=1))  # katal, catalitic activity
    # Non-SI but almost:
    new_mag('b', Magnitude(8.0, b=1))  # byte, note that B is Bel, as in dB
    
    # Funny Magnitudes with strange names
    # length
    new_mag("'", Magnitude(0.3048, m=1))  # feet
    new_mag('ft', Magnitude(0.3048, m=1))  # feet
    new_mag('inch', Magnitude(0.0254, m=1))  # inch
    new_mag('"', Magnitude(0.0254, m=1))  # inch
    new_mag('lightyear', Magnitude(2.99792458e8 * 365.25 * 86400, m=1))

    # volume
    new_mag('l', Magnitude(0.001, m=3))

    # time
    # year is tropical year, "the mean interval between vernal
    # equinoxes.  Differs from the sidereal year by 1 part in 26000
    # due to precession of the earth about its rotational axis
    # combined with precession of the perihelion of the earth's orbit"
    # (from units.dat).
    new_mag('year', Magnitude(31556925.974678401, s=1)) 
    new_mag('day', Magnitude(86400, s=1))
    new_mag('h', Magnitude(3600, s=1))
    new_mag('min', Magnitude(60, s=1))

    # Resolution
    new_mag('dpi', Magnitude(1.0 / 0.0254, m=-1))
    new_mag('lpi', Magnitude(1.0 / 0.0254, m=-1))

    # Velocity
    new_mag('ips', Magnitude(0.0254, m=1, s=-1))
    new_mag('c', Magnitude(2.99792458e8, m=1, s=-1))

    # Acceleration
    new_mag('gravity', Magnitude(9.80665, m=1, s=-2))


if not _mags:
    _init_mags()

if __name__ == '__main__':
    import doctest
    doctest.testmod()
