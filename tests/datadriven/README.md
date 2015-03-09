
The tests in this directory are data-driven tests, and use a custom 
extension to nose's test generators. The basic idea is that 
each test method that we define will be applied independently to 
each available data set. The data used to drive each of these tests
is found in the `tests/data/<data type>` directory. To add a new test,
simply add a new method that starts with `test` to a tester class.
To add new data, add the appropriate files into the `tests/data`
hierarachy.

The tests are automatically run using the nosetests. To run only 
the data driven tests (with verbose output) use
```
$ nosetests -v tests/datadriven
```
