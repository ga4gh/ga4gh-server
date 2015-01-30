# GA4GH Server Style Guide

The code follows the guidelines of `PEP 8
<http://legacy.python.org/dev/peps/pep-0008>`_ in most cases. The only notable
difference is the use of camel case over underscore delimited identifiers; this
is done for consistency with the GA4GH API. Code should be checked for compliance
using the `pep8 <https://pypi.python.org/pypi/pep8>`_ tool.

## Naming
Short variable names are OK when the scope is small (less than five lines) and the meaning is obvious. If the scope is longer or more complex, please use descriptive longer names.

## Imports
Imports should be structured into the following groups:

1. Any ```__future__``` imports
2. Any standard library imports
3. Any third party library imports
4. Any package level imports

## Tests
Tests should avoid use of the `assert` keyword and instead use the methods that `unittest` provides.
