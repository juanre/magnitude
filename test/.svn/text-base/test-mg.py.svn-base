#!/usr/bin/env python

"""
>>> from magnitude import mg, Magnitude, new_mag
>>> m = mg(1, '1/s')
>>> (m.val, m.unit) == (1.0, [0, -1, 0, 0, 0, 0, 0, 0, 0])
True
>>> lb = Magnitude(0.45359237, kg=1)
>>> new_mag('lb', lb)
>>> (lb.val - 0.45359237) < 1e-6
True
>>> lb.unit == [0, 0, 0, 1, 0, 0, 0, 0, 0]
True
>>> lb2 = mg(1, 'lb')
>>> lb2 == lb
True
>>> m = mg(10, 'm2 K / min')
>>> m.val - 10.0/60.0 < 1e-6
True
>>> m.unit
[2, -1, 1, 0, 0, 0, 0, 0, 0]
>>> m = mg(10, 'm/s') + (1, 'lightyear/year')
>>> n = mg(10, 'm/s') + mg(1, 'lightyear/year')
>>> o = (1, 'lightyear/year') + mg(10, 'm/s')
>>> p = (5, 'm/s') + mg(1, 'lightyear/year') + (5, 'm/s')
>>> m == n
True
>>> n == o
True
>>> p == m
True
>>> m.unit
[1, -1, 0, 0, 0, 0, 0, 0, 0]
>>> m = mg(10, 'km/h') - (1, 'km/h')
>>> n = (10, 'km/h') - mg(1, 'km/h')
>>> m.unit
[1, -1, 0, 0, 0, 0, 0, 0, 0]
>>> m.val == 10000.0 / 3600 - 1000.0 / 3600
True
>>> m == n
True
>>> m = mg(10, 'm/s') * (2, 's')
>>> n = mg(20, 'm')
>>> m == n
True
>>> (10, 'm/s') * mg(2, 's') == n
True
>>> mg(10, 'm/s') / (0.5, 'm') == mg(20, '1/s')
True
>>> (10, 'm/s') / mg(0.5, 'm') == mg(20, '1/s')
True
>>> 1 / mg(2, 's') == mg(0.5, '1/s')
True
>>> mg(5, 'kg/m2') % 3 == mg(2, 'kg/m2')
True
>>> mg(5, 'm/s') // mg(3, 'm') == mg(1, '1/s')
True
>>> mg(10, 'm/s') ** 2 == mg(100, 'm2/s2')
True
>>> mg(100, 'm2/s2') ** 0.5 == mg(10, 'm/s')
True
>>> (-mg(10, 'm/s')).val == -10.0
True
>>> (+mg(10, 'm/s')).val == 10.0
True
>>> abs(-mg(10, 'm/s')) == mg(10, 'm/s')
True
>>> int(mg(10.9, 'm/s')) == 10
True
>>> float(mg(10.9, 'm/s')) == 10.9
True
>>> long(mg(10.9, 'm/s')) == 10
True
>>> m = mg(10, 'm/s')
>>> m += mg(3, 'm/s')
>>> m.val
13.0
>>> m -= mg(43, 'm/s')
>>> m.val
-30.0
>>> m *= mg(2, 's')
>>> m == mg(-60, 'm')
True
>>> m = mg(10, 'm2/s')
>>> m /= (5, 'm')
>>> m == mg(2, 'm/s')
True
>>> m = mg(5, 'm/s')
>>> m //= mg(3, 'm')
>>> m == mg(1, '1/s')
True
>>> m = mg(5, 'm/s')
>>> m **= 2
>>> m == mg(25, 'm2/s2')
True
>>> mg(1, 'm/s') != mg(2, 'm/s')
True
>>> mg(1, 'm/s') < mg(2, 'm/s')
True
>>> mg(1, 'm/s') > mg(2, 'm/s')
False
>>> a = mg(1000, 'm/s')
>>> b = a.ounit('km/s')
>>> b.toval()
1.0
>>> m = mg(10, 'Gibit')
>>> m.val
10737418240.0
"""

### Tests I'm going to ignore, at least for the time being.
##     print mg(1, 'pl / [300]')
    
##     # mod
##     print (5, 'kg/m') % mg(3, 'kg/m2')
##     # divmod
##     print "divmod(10 m/s, 3 m) -> ", divmod(mg(10, 'm/s'), (3, 'm'))[0], \
##           ', ', divmod(mg(10, 'm/s'), (3, 'm'))[1]
##     # rdivmod
##     print "divmod(10 m/s, 3 m) -> ", divmod((10, 'm/s'), mg(3, 'm'))[0], \
##           ', ', divmod((10, 'm/s'), mg(3, 'm'))[1]
##     print "(10 m/s ** 2) % 3 -> ", pow(mg(10, 'm/s'), 2, 3)
##     print "10 m2/s2 ** 1/2 -> ", mg(10, 'm2/s2') ** 0.5
##     # imod
##     a = mg(10, 'm/s'); a %= (3, 'm/s')
##     print "a = mg(10, 'm/s'); a %= (3, 'm/s') -> ", a

if __name__ == '__main__':
    import doctest
    doctest.testmod()
