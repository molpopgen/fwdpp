# fwdpp testing suite

Unit and integration tests are found in the unit subdirectory of the source code repository.

The source code for unit tests are in unit/ and integration tests in integration/.

These tests make use of
[fixtures](http://www.boost.org/doc/libs/1_60_0/libs/test/doc/html/boost_test/tests_organization/fixtures/case.html) so
that a common set of objects can be re-used for similar tests.

##Dependencies

1. The [boost](http://boost.org) unit testing library is used by these tests. 

##Compiling the tests

~~~
make check
~~~

##Running the tests

~~~~
make check
~~~~

If you really want all the details, then execute this instead:

~~~
BOOST_TEST_LOG_LEVEL=all make check
~~~

The boost unit testing library will report any errors in any testing modules.

Note that some tests may intentionally cause errors.  When that it the case, a message stating that the error is intentional will appear on screen along with the error.

##Notes

* More tests are needed!
