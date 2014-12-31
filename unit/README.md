# Unit tests for fwdpp

##Dependencies

1. The [boost](http://boost.org) unit testing library is used by these tests. 

##Compiling the tests

```
make check
```

##Running the tests

```
sh runTests.sh
```

If you really want all the details, then execute this instead:

```
BOOST_TEST_LOG_LEVEL=all sh runTests.sh
```

The boost unit testing library will report any errors in any testing modules.

Note that some tests may intentionally cause errors.  When that it the case, a message stating that the error is intentional will appear on screen along with the error.

##Notes

* More tests are needed!
