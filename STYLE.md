# GA4GH Server Style Guide

The code follows the guidelines of [PEP 8]
(<http://legacy.python.org/dev/peps/pep-0008>) in most cases. The only notable
difference is the use of camel case over underscore delimited identifiers; this
is done for consistency with the GA4GH API. Code should be checked for compliance
using the [pep8](<https://pypi.python.org/pypi/pep8>) tool.

## Naming
Short variable names are OK when the scope is small (less than five lines)
and the meaning is obvious. If the scope is longer or more complex,
please use descriptive longer names.

Contrary to PEP8, acronyms should be camel cased, not capitalized.  For instance, HTTPAPIClient should be written as HttpApiClient.

## Imports
Imports should be structured into the following groups:

1. Any ```__future__``` imports
2. Any standard library imports
3. Any third party library imports
4. Any package level imports

The preferred approach is to import modules, and use this module
as a prefix to make it clear where a function or
class has come from. For example:
```python
# Good
import itertools
iterator = itertools.combinations(range(4), 2)

# Bad
from itertools import combinations
iterator = combinations(range(4), 2)
```

For deeper imports where it is not practical to type out the fully
qualified package name each time, it's useful to use the last 
part of the package name as a shortcut. For example:
```python
# Good
import ga4gh.protocol as protocol
variant = protocol.Variant()

# Not so good 
import ga4gh
variant = ga4gh.protocol.Variant()

# Bad
from ga4gh.protocol import Variant
variant = Variant()
```

## Tests
Tests should avoid use of the `assert` keyword and instead use the methods that `unittest` provides.
