To run the unit tests, run `./fetch_gtest.sh && make test`, then run each of the
resulting test executables that get produced in the `bin/` directory.

To add a new file full of tests "foo_test.cpp", add bin/foo_test to the
Makefile's `test:` target, then add foo_test.cpp to the `test/` directory.
Note that the name must end in _test for `make` to know how to handle it.
