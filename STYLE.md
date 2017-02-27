# GA4GH Server Style Guide

The code follows the guidelines of [PEP 8]
(<http://legacy.python.org/dev/peps/pep-0008>) in most cases. The only notable
difference is the use of camel case over underscore delimited identifiers; this
is done for consistency with the GA4GH API. Code should be checked for compliance
using the [pep8](<https://pypi.python.org/pypi/pep8>) tool.  In practice, we are
a bit inconsistent about this since switching over to protobuf, which uses
snake_cased names.

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
4. Any package level imports from the current package
5. Any imports from other ga4gh packages

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
Unit tests should have descriptive 
method names and should aim to be short, simple and isolated. Where complex 
testing is required, we should delegate the verification to another 
method. For example, 

```python
class TestGenotypeConversion(unittest.TestCase):
    """ 
    Tests the conversion of VCF genotypes to GA4GH values.
    """
    def verifyConversion(self, vcfGenotype, vcfPhaseset, 
                         ga4ghGenotype, ga4ghPhaseset):
        """
        Verifies that the convertGenotype function properly converts the 
        specified VCF genotype and phaseset values into the specified 
        GA4GH genotype and phaseset values.
        """
        self.assertEqual(
            (ga4ghGenotype, ga4ghPhaseset),
            variants.convertVcfGenotype(vcfGenotype, vcfPhaseset))

    def testUnphasedNoCall(self):
        self.verifyConversion("./.", "0", [-1], None)

    def testUnphasedSecondHalfCall(self):
        self.verifyConversion("./0", "25", [-1], None)

    # etc.
```

Tests should avoid use of the `assert` keyword and instead use the
methods that `unittest` provides. 

## Docstring comments

All public methods and functions should have docstring comments, and should 
follow the conventions of [PEP 257](https://www.python.org/dev/peps/pep-0257/).

The exception to this rule is for unit tests, which should have descriptive
method names, and no docstring comments. (The reason for this is that 
nosetests prints out the docstring comment rather than the methodname 
when running `nosetests -v`, which seems less useful.) See above for an 
example of the recommended structure for unit tests. 


