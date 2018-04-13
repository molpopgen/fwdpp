# Tutorial 1: Policies in fwdpp

[TOC]

\section TutIntro Introduction

Please read @ref md_md_datatypes and @ref md_md_details before reading this document.

This document is intended to be an in-depth discussion of policies and their role in implementing forward-time population genetic simulations using the C++ template library __fwdpp__.  We will first describe what policies are using standard C++ examples, and then we will get into the harder stuff.

An understanding of C++ fundamentals, including containers, their iterators, and how they relate to the standard algorithms, is assumed knowledge here.

\section TutGeneral Policies in C++

The behavior of a function may be modified by a policy.  For example, [std::sort] defaults to sorting in ascending
order:

~~~{.cpp}
std::vector<int> x{5,4,3,2,1};
std::sort(x.begin(),x.end()); //x will now contain 1,2,3,4,5
~~~

This default behavior may be modified by passing a custom function ("policy") as an additional argument.  In C++, this
function can take on several different forms, including function objects, lambda expression, and instances of
std::function. Here is an _incomplete_ list of the various ways you could sort in ascending order:

~~~{.cpp}
#include <algorithm>
#include <vector>
#include <functional>

using namespace std;

int main(int argc, char ** argv) {
	vector<int> x{5,4,3,2,1};

	sort(x.begin(),x.end()); //x is now 1,2,3,4,5

	//sort x in ascending order via std::greater...
	sort(x.begin(),x.end(),bind(greater<int>(),placeholders::_1,placeholders::_2));

	//...or a lambda function...
	sort(x.begin(),x.end(),[](int a,int b) { return a>b; });

	//...or a std::function wrapping a lambda...
	std::function<bool(int,int)> sorter = [](int a,int b){return a>b;};
	sort(x.begin(),x.end(),sorter);

	//...or a std::function wrapping bind/greater...
	sorter = std::bind(greater<int>(),placeholders::_1,placeholders::_2);
	sort(x.begin(),x.end(),sorter);
}
~~~

\section fwdppPolicies Policies in fwdpp

Just like the STL, custom functions modify the behavior of simulations.  These policies affect how mutation,
recombination, and selection all work.  Additionally, policies may affect how mutations are removed from gametes and/or
the entire population at the end of each generation.  For example, consider a neutral mutation that fixed in the most
recent generation.  There is really no point in keeping a key to it in every gamete.  Rather, we may prefer to remove it
from each gamete and store it in a vector of fixations.  Doing so requires a policy passed to fwdpp::sample_diploid
telling it to remove neutral mutations.  In the case of selected mutations, we may or may not want to remove such fixed
variants from gametes, depending on the context of our simulation (standard pop-gen where fixed selective mutations no
longer affect relative fitness	vs quantitative trait simulations where they do contribute to trait values, and thus
affect distance to optimium, etc.).

\subsection TutMut Mutation policies

A mutation model has one of the following function signatures. Either 

~~~{.cpp}
std::size_t mutmodel(queue_t &,mcont_t &);
~~~

or

~~~{.cpp}
std::size_t mutmodel(gamete_t &,mcont_t &, queue_t &);
~~~

The former is used when the current state of a gamete doesn't matter.  The infinitely-many sites model is an example.
The latter is used whn the current state of a gamete does matter.  A finite-sites model would be an example.
Internally, the library detects which type is being used.  The return value is the index in the mutation container where
the new mutation was placed.  This location may have been generated via recycling, and the library provides
fwdpp::fwdpp_internal::recycle_mutation_helper to facilitate this operation.

Additional model parameters are possible via the usual mechanisms: class members, std::bind, etc.

A valid mutation policy passes a static assertion involving fwdpp::traits::valid_mutation_model at compile time. See
type_traitsTest.cc for an example of this assertion.

For an example of how to compose a mutation policy, see @ref custom_mutation.cc.

\subsection TutRec Recombination policies

The function signature of a recombination policy must be equivalent to fwdpp::traits::recmodel_t.  If so, then a static
assertion involving fwdpp::traits::valid_rec_model will pass at compile time. See type_traitsTest.cc for an example of
this assertion. 

Additional model parameters are possible via the usual mechanisms: class members, std::bind, etc.

* fwdpp::poisson_xover

\subsection TutFitness Genetic value policies

* fwdpp::site_dependent_genetic_value
* fwdpp::additive_diploid
* fwdpp::multiplicative_diploid

\subsection TutRemoval "Removal policies"

The default is std::true_type, which means that all fixations are removed from gametes.

* fwdpp::remove_nothing
* fwdpp::remove_neutral
