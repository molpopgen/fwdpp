/*
Essentially a re-implementation of diploid_ind.cc that is capable of being run
in python.

This implementation is in terms of pybind11: https://github.com/pybind/pybind11

In addition to the usual fwdpp depdencies, we need pybind11 installed

To compile:
g++ -fPIC -Wall -W -O3 -I. `python-config --includes` -std=c++11 -c
fwdpp_pybind11.cc
g++ -std=c++11 -shared -o fwdpp_pybind11.so fwdpp_pybind11.o -lpython -lgsl
-lgslcblas

To run:
python test_pybind11.py
*/

#include <python_common.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Calculate the site-frequency spectrum for a sample
std::vector<unsigned>
sfs(GSLrng &rng, const poptype &pop, const unsigned &nsam)
{
    std::map<double, unsigned> mutfreqs;
    unsigned twoN = 2 * pop.N;

    for (unsigned i = 0; i < nsam; ++i)
        {
            // pick a random chrom (w/replacement...)
            unsigned chrom
                = unsigned(gsl_ran_flat(rng.get(), 0., double(twoN)));
            // get reference to that chrom from the individual
            auto &gamete = (chrom % 2 == 0.)
                               ? pop.gametes[pop.diploids[chrom / 2].first]
                               : pop.gametes[pop.diploids[chrom / 2].second];
            // In this example, there are only neutral mutations, so that's
            // what we'll iterate over
            for (auto &m : gamete.mutations)
                {
                    auto pos = pop.mutations[m].pos;
                    auto pos_itr = mutfreqs.find(pos);
                    if (pos_itr == mutfreqs.end())
                        {
                            mutfreqs.insert(std::make_pair(pos, 1));
                        }
                    else
                        {
                            pos_itr->second++;
                        }
                }
        }
    // Now, fill in the SFS, omitting positions that are fixed in the sample
    std::vector<unsigned> __rv(nsam - 1, 0u);
    for (const auto &__x : mutfreqs)
        {
            if (__x.second < nsam)
                __rv[__x.second - 1]++;
        }
    return __rv;
}

// Now, we can expose the stuff to python
PYBIND11_PLUGIN(fwdpp_pybind11)
{
    pybind11::module m("fwdpp_pybind11",
                       "example of wrapping fwdpp in python using pybind11");

    /*
      Make the mutation types from fwdpp visible to Python.

      For a "real" extension, we'd do this more rigorously,
      and first expose KTfwd::mutation_base,
      and then tell pybind11 that KTfwd::mutation inherits from that.
    */
    pybind11::class_<KTfwd::mutation>(m, "mutation")
        .def_readonly("pos", &KTfwd::mutation::pos)
        .def_readonly("s", &KTfwd::mutation::s)
        .def_readonly("h", &KTfwd::mutation::h)
        /*
          This is where pybind is really cool.
          We can write a lambda here to
          get a representation of this type as a
          Python dictionary.
        */
        .def("as_dict", [](const KTfwd::mutation &m) noexcept {
            using obj = pybind11::object;
            pybind11::dict rv;
            rv[obj(pybind11::cast("pos"))] = obj(pybind11::cast(m.pos));
            rv[obj(pybind11::cast("s"))] = obj(pybind11::cast(m.s));
            rv[obj(pybind11::cast("h"))] = obj(pybind11::cast(m.h));
            return rv;
        });
    ;

    // Expose gamete type
    pybind11::class_<KTfwd::gamete>(m, "gamete")
        .def_readonly("n", &KTfwd::gamete::n)
        .def_readonly("mutations", &KTfwd::gamete::mutations)
        .def_readonly("smutations", &KTfwd::gamete::smutations)
        .def("as_dict",
             // This lambda shows that
             // we can return dicts
             // with a mix of scalars + containers
             [](const KTfwd::gamete &g) noexcept {
                 using obj = pybind11::object;
                 pybind11::dict rv;
                 rv[obj(pybind11::cast("n"))] = obj(pybind11::cast(g.n));
                 rv[obj(pybind11::cast("mutations"))]
                     = obj(pybind11::cast(g.mutations));
                 rv[obj(pybind11::cast("smutations"))]
                     = obj(pybind11::cast(g.smutations));
                 return rv;
             });
    ;

	using poptype_base = poptype::popbase_t;
	
	pybind11::class_<poptype_base>(m,"poptype_base")
	;

    // Expose the type based on fwdpp's "sugar" layer
    pybind11::class_<poptype,poptype_base>(m, "poptype")
        .def(pybind11::init<unsigned>())
        .def("clear", &poptype::clear)
        /*
          The next three show how to
          use pybind's auto-conversion of std::vector
          to Python list types.  This is a big improvement
          over boost.python!
        */
        .def_readonly("mutations", &poptype::mutations)
        .def_readonly("mcounts", &poptype::mcounts)
        .def_readonly("fixations", &poptype::fixations)
        .def_readonly("diploids", &poptype::diploids)
        .def_readonly("gametes", &poptype::gametes);

    // Expose the GSL wrapper
    pybind11::class_<GSLrng>(m, "GSLrng").def(pybind11::init<unsigned>());
    // Expose the function to run the model
    // Note the use of C++11 string literals
    // to make a nice docstring that Sphinx can parse out.
    m.def("evolve", &evolve, R"delimiter(Evolve a population

:param rng: A GSLrng
:param N: Diploid population size
:param generations: Number of generations to simulation
:param mu: Mutation rate (per gamete, per generation)
:param recrate: Recombination rate (per diploid, per generation))delimiter");
    // And one to get the sfs of a sample
    m.def("sfs", &sfs);

    return m.ptr();
}
