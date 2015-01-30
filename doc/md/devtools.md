# Tutorial 4: Starting projects based on fwdpp

## Getting the right tools on OS X

OS X users are at a slight disadvantage here, as the Xcode tool suite no longer includes the GNU tools setup.  However, you can use [homebrew](http://brew.sh) to get these tools.  Please see that website for complete instructions, but at the time of this writing, the following will work:

1. Make sure that you have Xcode installed.  You cannot in stall brew without it.
2. Install brew:

~~~{.sh}
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
~~~

3. Install several packages that you'll need:

~~~{.sh}
brew install autoconf automake libtool wget
~~~

4. Optionally, install some other handy tools:

~~~{.sh}
brew install doxygen
~~~
