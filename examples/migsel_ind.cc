/*! \include migsel_ind.cc
  Two constant-size populations with selection and inbreeding.

  The selection coefficient, s, is treated as -s in deme 2, just for fun.

  Writes the metapop + an "ms"-type sample in binary format to an output file.
*/
#include <config.h>
#include <numeric>
#include <cmath>
#include <functional>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <fwdpp/diploid.hh>

#ifdef HAVE_LIBSEQUENCE
#include <Sequence/SimData.hpp>
#include <Sequence/SimDataIO.hpp> //for writing & reading SimData objects in binary format
#include <Sequence/FST.hpp>
#endif
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/serialization.hpp>
// the type of mutation
using mtype = KTfwd::mutation;
#define METAPOP_SIM
#include <common_ind.hpp>

using poptype = metapop_t;
using mcont = poptype::mcont_t;
using gtype = poptype::gamete_t;
using gcont = poptype::gcont_t;

using namespace std;
using namespace KTfwd;
#ifdef HAVE_LIBSEQUENCE
using namespace Sequence;
#endif

size_t
migpop(const size_t &source_pop, const gsl_rng *r, const double &mig_prob)
{
    if (gsl_rng_uniform(r) < mig_prob)
        {
            return !source_pop;
        }
    return source_pop;
}

#ifdef HAVE_LIBSEQUENCE
SimData merge(const std::vector<std::pair<double, std::string>> &sample1,
              const std::vector<std::pair<double, std::string>> &sample2,
              const unsigned &nsam);
#endif

// fitness model details -- s will be treated as -s in population 2
struct multiplicative_diploid_minus
{
    typedef double result_type;
    inline double
    operator()(const poptype::diploid_t &dip, const poptype::gcont_t &gametes,
               const poptype::mcont_t &mutations,
               const double scaling = 1.) const
    {
        using mut_t = poptype::mcont_t::value_type;
        return site_dependent_fitness()(
            gametes[dip.first], gametes[dip.second], mutations,
            [&](double &fitness, const mut_t &mut) {
                fitness *= (1. - scaling * mut.s);
            },
            [](double &fitness, const mut_t &mut) {
                fitness *= (1. - mut.h * mut.s);
            },
            1.);
    }
};

int
main(int argc, char **argv)
{
    if (argc != 14)
        {
            std::cerr << "Too few arguments.\n"
                      << "Usage: migsel_ind N theta_neutral theta_deleterious "
                         "rho M s h f1 f2 ngens n outfilename seed\n";
            exit(0);
        }
    int argn = 1;
    const unsigned N = atoi(argv[argn++]);
    const double theta_neut = atof(argv[argn++]);
    const double theta_del = atof(argv[argn++]);
    const double rho = atof(argv[argn++]);
    const double M = atof(argv[argn++]);
    const double s = atof(argv[argn++]);
    const double h = atof(argv[argn++]);
    const double f1 = atof(argv[argn++]);
    const double f2 = atof(argv[argn++]);
    const unsigned ngens = atoi(argv[argn++]);
    const unsigned n = atoi(argv[argn++]);
    const char *outfilename = argv[argn++];
    const unsigned seed = atoi(argv[argn++]);

    const double mu_neutral = theta_neut / double(4 * N);
    const double mu_del = theta_del / double(4 * N);
    const double littler = rho / double(4 * N);
    const double m = M / double(4 * N);

    GSLrng r(seed);

    unsigned Ns[2] = { N, N };
    double fs[2] = { f1, f2 };
    poptype pop({ N, N });
    pop.mutations.reserve(
        2 * size_t(std::ceil(std::log(2 * N) * (theta_neut + theta_del)
                             + 0.667 * (theta_neut + theta_del))));
    // create a vector of fitness functions for each population
    std::vector<std::function<double(const poptype::diploid_t &,
                                     const poptype::gcont_t &,
                                     const poptype::mcont_t &)>>
        vbf;
    vbf.push_back(std::bind(multiplicative_diploid(), std::placeholders::_1,
                            std::placeholders::_2, std::placeholders::_3, 2.));
    vbf.push_back(std::bind(multiplicative_diploid_minus(),
                            std::placeholders::_1, std::placeholders::_2,
                            std::placeholders::_3, 2.));

    for (unsigned generation = 0; generation < ngens; ++generation)
        {
            std::vector<double> wbars = sample_diploid(
                r.get(), pop.gametes, pop.diploids, pop.mutations, pop.mcounts,
                &Ns[0], mu_neutral + mu_del,
                std::bind(KTfwd::infsites(), std::placeholders::_1,
                          std::placeholders::_2, r.get(),
                          std::ref(pop.mut_lookup), mu_neutral, mu_del,
                          [&r]() { return gsl_rng_uniform(r.get()); },
                          [&s]() { return s; }, [&h]() { return h; }),
                std::bind(KTfwd::poisson_xover(), r.get(), littler, 0., 1.,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                vbf, std::bind(migpop, std::placeholders::_1, r.get(), m),
                pop.neutral, pop.selected, &fs[0]);
            // 4*N b/c it needs to be fixed in the metapopulation
            update_mutations(pop.mutations, pop.fixations, pop.fixation_times,
                             pop.mut_lookup, pop.mcounts, generation, 4 * N);
        }

    std::pair<std::vector<std::pair<double, std::string>>,
              std::vector<std::pair<double, std::string>>>
        spop1 = ms_sample_separate(r.get(), pop.mutations, pop.gametes,
                                   pop.diploids[0], n);

    std::pair<std::vector<std::pair<double, std::string>>,
              std::vector<std::pair<double, std::string>>>
        spop2 = ms_sample_separate(r.get(), pop.mutations, pop.gametes,
                                   pop.diploids[1], n);

    std::ofstream outstream(outfilename);

    // Write the metapop in binary format to outstream
    KTfwd::write_binary_metapop(pop.gametes, pop.mutations, pop.diploids,
                                std::bind(KTfwd::mutation_writer(),
                                          std::placeholders::_1,
                                          std::placeholders::_2),
                                outstream);

    // Write the "ms" blocks
#ifdef HAVE_LIBSEQUENCE
    SimData neutral = merge(spop1.first, spop2.first, n);
    SimData selected = merge(spop1.second, spop2.second, n);
    Sequence::SimData neutral2, selected2;
    Sequence::write_SimData_binary(outstream, neutral);
    Sequence::write_SimData_binary(outstream, selected);
#endif
    outstream.close();

    poptype::gcont_t metapop2;
    poptype::vdipvector_t diploids2;
    poptype::mcont_t mutations2;

    ifstream in(outfilename);

    KTfwd::read_binary_metapop(
        metapop2, mutations2, diploids2,
        std::bind(KTfwd::mutation_reader<mtype>(), std::placeholders::_1), in);

    assert(metapop2.size() == pop.gametes.size());
    assert(mutations2.size() == pop.mutations.size());
    assert(diploids2.size() == pop.diploids.size());

#ifdef HAVE_LIBSEQUENCE
    neutral2 = Sequence::read_SimData_binary(in);
    selected2 = Sequence::read_SimData_binary(in);
    std::cerr << (neutral == neutral2) << ' ' << (selected == selected2)
              << '\n';
#endif
    in.close();

    /*
      At this point, you could go through each deme and each diploid and make
      sure that all is cool.  However, if we weren't reading and
      writing the metapop correctly, there's no way we'd be able
      to write and then read in the ms blocks correctly, as we'd have
      run into some binary gibberish along the way.
    */

    // For fun, we'll calculate some basic pop subdivision stats
#ifdef HAVE_LIBSEQUENCE
    unsigned config[2] = { n, n };
    if (!neutral.empty())
        {
            Sequence::FST fst_neut(&neutral, 2, config);
            std::pair<std::set<double>, std::set<double>> pneut
                = fst_neut.Private(0, 1);
            std::cout << fst_neut.HSM() << '\t' << fst_neut.shared(0, 1).size()
                      << '\t' << pneut.first.size() << '\t'
                      << pneut.second.size() << '\t';
        }
    else
        {
            std::cout << "NA\t0\t0\t0\t0\t";
        }
    if (!selected.empty())
        {
            Sequence::FST fst_sel(&selected, 2, config);
            std::pair<std::set<double>, std::set<double>> psel
                = fst_sel.Private(0, 1);
            std::cout << fst_sel.HSM() << '\t' << fst_sel.shared(0, 1).size()
                      << '\t' << psel.first.size() << '\t'
                      << psel.second.size() << '\n';
        }
    else
        {
            std::cout << "NA\t0\t0\t0\t0\n";
        }
#endif
}

#ifdef HAVE_LIBSEQUENCE
SimData
merge(const std::vector<std::pair<double, std::string>> &sample1,
      const std::vector<std::pair<double, std::string>> &sample2,
      const unsigned &nsam)
{
    std::map<double, std::string> temp;

    for (unsigned i = 0; i < sample1.size(); ++i)
        {
            temp[sample1[i].first]
                = std::string(sample1[i].second + std::string(nsam, '0'));
        }

    for (unsigned i = 0; i < sample2.size(); ++i)
        {
            std::map<double, std::string>::iterator itr
                = temp.find(sample2[i].first);
            if (itr == temp.end())
                {
                    temp[sample2[i].first] = std::string(std::string(nsam, '0')
                                                         + sample2[i].second);
                }
            else
                {
                    std::copy(sample2[i].second.begin(),
                              sample2[i].second.end(),
                              itr->second.begin() + nsam);
                }
        }
    std::vector<std::pair<double, std::string>> rv(temp.begin(), temp.end());
    std::sort(rv.begin(), rv.end(), [](std::pair<double, std::string> lhs,
                                       std::pair<double, std::string> rhs) {
        return lhs.first < rhs.first;
    });
    return SimData(rv.begin(), rv.end());
}
#endif
