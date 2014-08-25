bin_PROGRAMS=FWDPPPACKAGE

FWDPPPACKAGE_SOURCES=FWDPPPACKAGE.cc

AM_CPPFLAGS=-W -Wall -I.

if DEBUG
else !DEBUG
AM_CPPFLAGS+=-DNDEBUG
endif
LDADD=

if BVEC_PRESENT
if BLIST_PRESENT
if BPOOL_PRESENT
if BHASH_PRESENT
if BSET_PRESENT
#On some systems (or at least recent boost versions, the boost_system is needed to use the boost containers
endif
endif
endif
endif
else
#Boost is missing, so we gotta use the standard headers                                                                   
AM_CPPFLAGS+=-DUSE_STANDARD_CONTAINERS
endif
