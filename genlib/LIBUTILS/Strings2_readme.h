STRINGS, version 1.2
Mart Rentmeester, Mart.Rentmeester@nn-online.org
http://nn-online.org/code/



STRINGS

This STRINGS module is very much inspired by the 'iso_varying_strings' module
and the sample implementation by L. Schonfelder. Most of the functions (except
I/O functions) mentioned in the ISO/IEC 1539-2:1999 extension to the fortran 95
standard are implemented one way or the other. Some other useful routines have
been added to the module, like sorting arrays of character variables or
strings, and transforming to uppercase/lowercase.

The main difference with the iso_varying_strings module however is that the
functions in this strings module return a value of fortran's native character
type instead of 'string' type. Although that usually may lead to some
performance penalty due to extra copying, I see two main advantages to this
approach:

  * It is more useful. The character variable in fortran is very powerful. And
    the function results can be used directly in output statements.
  * It is possible to write the module (if I have done everything correctly)
    without memory leaks in standard fortran 95. New allocation only takes
    place during assignment of local string variables. It is always up to the
    user to deallocate used local string variables using the 'unstring'
    subroutine!!!!!! (It is possible to compile the routine to use allocatable
    components as described in technical report 15581 instead of using
    pointers. In that case there may no longer be a need to use the unstring
    routine.)

In addition there is the possibility to play with the length of the allocated
space that holds the string as if it were a kind of buffer. Several subroutines
to manipulate the string and the allocated space are provided.

Alas no documentation yet. With a bit of luck near the end of 2003.

The code is provided as is, you can use it as you wish. But I would appreciate
it if you'd let me know if you find a particularly useful purpose for it - I
might make some kind of gallery. Feedback, suggestions and contributions to
improve this software are more than welcome; if you have any questions,
remarks, bug fixes, donations, postcards, or job offers, please contact me:
Mart.Rentmeester@nn-online.org

