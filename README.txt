Copyright (C) 2006-2015 Juan Reyero (http://juanreyero.com).

Licensed under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
either express or implied.  See the License for the specific
language governing permissions and limitations under the
License.

Home page: http://juanreyero.com/open/magnitude/

About:

Magnitude is a library for computing with numbers with
units::

  >>> print "10 m/s ** 2 ->", mg(10, 'm/s') ** 2
  10 m/s ** 2 -> 100.0000 m2 / s2
  >>> print (mg(10, 'm') * 2 / (10, 'm/s2')).sqrt()
  1.4142 s
  >>> tsq = mg(10, 'm') * 2 / (10, 'm/s2')
  >>> print tsq ** 0.5
  1.4142 s
  >>> print mg(1, "lightyear") / mg(1, "c")
  >>> 31557600.0000 s
  >>> y = mg(1, "lightyear") / (1, "c")
  >>> y.ounit("year")
  <magnitude.magnitude instance at 0x81440>
  >>> print y
  1.0000 year
  >>> yd = y.ounit('day')
  >>> print yd
  365.2500 day
  >>> output_precision(0)
  0
  >>> print yd
  365 day
