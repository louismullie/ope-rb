##ope-rb

This gem implements order-preserving symmetric encryption, as described in [Boldyreva, 2009](http://www.cc.gatech.edu/~aboldyre/papers/bclo.pdf).

###Usage

```ruby

require 'ope'

cipher = OPE::Cipher.new(key)
a = cipher.encrypt(456789)
b = cipher.encrypt(891234)

puts b > a     # =>  true
```

###Credits

[CryptDB](http://g.csail.mit.edu/cryptdb/) provided a reference implementation for the Boldyreva paper. The code from [Caesar](https://github.com/Bren2010/caesar) was helpful in understanding how the scheme worked from a high-level perspective. The hypergeometric distribution code was adapted from a [Fortran implementation](http://calgo.acm.org/) from the Association for Computing Machinery and a [C implementation](http://ics.hutton.ac.uk/svn/topali-v2/trunk/binaries/src/barce/mathlib.cc) used by the R project.

###License

This software is released under the GNU Affero General Public License. If this license does not suit your needs, please contact me.