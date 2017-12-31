#ifndef __FWDPP_SUGAR_SINGLEPOP_SERIALIZATON_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_SERIALIZATON_HPP__

#include <ios>
#include <iosfwd>
#include <sstream>
#include <algorithm>
#include <type_traits>
#include <utility>
#include <functional>
#include <fwdpp/IO.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/generalmut.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

namespace fwdpp
{
    namespace detail
    {
        template <typename streamtype, typename poptype>
        inline void
        serialize_population_details(streamtype &buffer, const poptype &pop,
                                     sugar::SINGLEPOP_TAG)
        {
            io::scalar_writer writer;
            writer(buffer, &pop.N);
            write_binary_pop(pop.gametes, pop.mutations, pop.diploids, buffer);
            // Step 2: output fixations
            fwdpp::io::write_mutations(pop.fixations, buffer);
            if (!pop.fixations.empty())
                {
                    // Step 3:the fixation times
                    writer(buffer, &pop.fixation_times[0],
                           pop.fixations.size());
                }
        }

        template <typename streamtype, typename poptype>
        inline void
        serialize_population_details(streamtype &buffer, const poptype &pop,
                                     sugar::MULTILOCPOP_TAG)
        {
            io::scalar_writer writer;
            writer(buffer, &pop.N);
            write_binary_pop_mloc(pop.gametes, pop.mutations, pop.diploids,
                                  buffer);
            // Step 2: output fixations
            fwdpp::io::write_mutations(pop.fixations, buffer);
            if (!pop.fixations.empty())
                {
                    // Step 3:the fixation times
                    writer(buffer, &pop.fixation_times[0],
                           pop.fixations.size());
                }
            std::size_t temp = std::size_t(pop.locus_boundaries.size());
            writer(buffer, &temp);
            for (std::size_t i = 0; i < temp; ++i)
                {
                    writer(buffer, &pop.locus_boundaries[i].first);
                    writer(buffer, &pop.locus_boundaries[i].second);
                }
        }

        template <typename streamtype, typename poptype>
        inline void
        serialize_population_details(streamtype &buffer, const poptype &pop,
                                     sugar::METAPOP_TAG)
        {
            io::scalar_writer writer;
            std::size_t npops = std::size_t(pop.Ns.size());
            writer(buffer, &npops);
            writer(buffer, &pop.Ns[0], npops);
            write_binary_metapop(pop.gametes, pop.mutations, pop.diploids,
                                 buffer);
            // Step 2: output fixations
            fwdpp::io::write_mutations(pop.fixations, buffer);
            if (!pop.fixations.empty())
                {
                    // Step 3:the fixation times
                    writer(buffer, &pop.fixation_times[0],
                           pop.fixations.size());
                }
        }

        template <typename streamtype, typename poptype>
        inline void
        deserialize_population_details(poptype &pop, streamtype &buffer,
                                       sugar::SINGLEPOP_TAG)
        {
            pop.clear();
            io::scalar_reader reader;
            // Step 0: read N
            reader(buffer, &pop.N);
            fwdpp::read_binary_pop(pop.gametes, pop.mutations, pop.diploids,
                                   buffer);

            // update the mutation counts
            fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                            pop.mcounts);
            fwdpp::io::read_mutations(pop.fixations, buffer);
            if (!pop.fixations.empty())
                {
                    pop.fixation_times.resize(pop.fixations.size());
                    reader(buffer, &pop.fixation_times[0],
                           pop.fixations.size());
                }

            // Finally, fill the lookup table:
            for (unsigned i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (pop.mcounts[i])
                        pop.mut_lookup.insert(pop.mutations[i].pos);
                }
        }

        template <typename streamtype, typename poptype>
        inline void
        deserialize_population_details(poptype &pop, streamtype &buffer,
                                       sugar::MULTILOCPOP_TAG)
        {
            pop.clear();
            io::scalar_reader reader;
            // Step 0: read N
            reader(buffer, &pop.N);
            fwdpp::read_binary_pop_mloc(pop.gametes, pop.mutations,
                                        pop.diploids, buffer);
            // update the mutation counts
            fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                            pop.mcounts);
            std::size_t temp;
            fwdpp::io::read_mutations(pop.fixations, buffer);
            if (!pop.fixations.empty())
                {
                    pop.fixation_times.resize(pop.fixations.size());
                    reader(buffer, &pop.fixation_times[0],
                           pop.fixations.size());
                }
            reader(buffer, &temp);
            if (temp)
                {
                    double x[2];
                    for (std::size_t i = 0; i < temp; ++i)
                        {
                            reader(buffer, &x[0], 2);
                            pop.locus_boundaries.emplace_back(x[0], x[1]);
                        }
                }

            // Finally, fill the lookup table:
            for (unsigned i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (pop.mcounts[i])
                        pop.mut_lookup.insert(pop.mutations[i].pos);
                }
        }

        template <typename streamtype, typename poptype>
        inline void
        deserialize_population_details(poptype &pop, streamtype &buffer,
                                       sugar::METAPOP_TAG)
        {
            pop.clear();
            io::scalar_reader reader;
            // Step 0: read N
            std::size_t numNs;
            reader(buffer, &numNs);
            pop.Ns.resize(numNs);
            reader(buffer, &pop.Ns[0], numNs);
            // Step 1: write the mutations, diploids, gametes to the stream
            fwdpp::read_binary_metapop(pop.gametes, pop.mutations,
                                       pop.diploids, buffer);
            // update the mutation counts
            fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                            pop.mcounts);
            fwdpp::io::read_mutations(pop.fixations, buffer);
            if (!pop.fixations.empty())
                {
                    pop.fixation_times.resize(pop.fixations.size());
                    reader(buffer, &pop.fixation_times[0],
                           pop.fixations.size());
                }

            // Finally, fill the lookup table:
            for (unsigned i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (pop.mcounts[i])
                        pop.mut_lookup.insert(pop.mutations[i].pos);
                }
        }
    }

    template <typename streamtype, typename poptype>
    inline void
    serialize_population(streamtype &buffer, const poptype &pop)
	/// Write a population in binary format to a stream.
	///
	/// \param buffer A model of std::ostream.
	/// \param pop A population.
    {
        detail::serialize_population_details(buffer, pop,
                                             typename poptype::popmodel_t());
    }

    template <typename streamtype, typename poptype>
    inline void
    deserialize_population(poptype &pop, streamtype &buffer)
	/// Read a binary representation of a population from a stream.
	/// 
	/// \param pop A population to fill from \a buffer
	/// \param buffer A model of std::istream.
    {
        detail::deserialize_population_details(pop, buffer,
                                               typename poptype::popmodel_t());
    }
}

#endif
