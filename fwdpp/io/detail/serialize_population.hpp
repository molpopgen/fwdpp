#ifndef FWDPP_IO_SERIALIZE_POPULATION_DETAIL_HPP__
#define FWDPP_IO_SERIALIZE_POPULATION_DETAIL_HPP__

#include <fwdpp/io/scalar_serialization.hpp>
#include <fwdpp/io/mutation.hpp>
#include <fwdpp/io/gamete.hpp>
#include <fwdpp/io/diploid.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

namespace fwdpp
{
    namespace io
    {
        namespace detail
        {
            template <typename streamtype, typename poptype>
            inline void
            serialize_population_details(streamtype &buffer,
                                         const poptype &pop,
                                         sugar::SINGLELOC_TAG)
            {
                io::scalar_writer writer;
                writer(buffer, &pop.N);
                io::write_mutations(buffer, pop.mutations);
                io::write_gametes(buffer, pop.gametes);
                io::write_diploids(buffer, pop.diploids);
                // Step 2: output fixations
                fwdpp::io::write_mutations(buffer, pop.fixations);
                if (!pop.fixations.empty())
                    {
                        // Step 3:the fixation times
                        writer(buffer, &pop.fixation_times[0],
                               pop.fixations.size());
                    }
            }

            template <typename streamtype, typename poptype>
            inline void
            serialize_population_details(streamtype &buffer,
                                         const poptype &pop,
                                         sugar::MULTILOC_TAG)
            {
                io::scalar_writer writer;
                writer(buffer, &pop.N);
                unsigned nloci = unsigned(pop.diploids[0].size());
                writer(buffer, &nloci);
                // write mutations
                io::write_mutations(buffer, pop.mutations);
                io::write_gametes(buffer, pop.gametes);
                unsigned ndips = unsigned(pop.diploids.size());
                writer(buffer, &ndips);
                io::serialize_diploid<
                    typename poptype::dipvector_t::value_type::value_type>
                    dipwriter;
                for (const auto &dip : pop.diploids)
                    {
                        for (const auto &genotype : dip)
                            {
                                dipwriter(buffer, genotype);
                            }
                    }
                // Step 2: output fixations
                fwdpp::io::write_mutations(buffer, pop.fixations);
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
            deserialize_population_details(poptype &pop, streamtype &buffer,
                                           sugar::SINGLELOC_TAG)
            {
                pop.clear();
                io::scalar_reader reader;
                // Step 0: read N
                reader(buffer, &pop.N);
                io::read_mutations(buffer, pop.mutations);
                io::read_gametes(buffer, pop.gametes);
                io::read_diploids(buffer, pop.diploids);

                // update the mutation counts
                fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                                pop.mcounts);
                fwdpp::io::read_mutations(buffer, pop.fixations);
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
                            pop.mut_lookup.emplace(pop.mutations[i].pos,
                                                   static_cast<uint_t>(i));
                    }
            }

            template <typename streamtype, typename poptype>
            inline void
            deserialize_population_details(poptype &pop, streamtype &buffer,
                                           sugar::MULTILOC_TAG)
            {
                pop.clear();
                io::scalar_reader reader;
                // Step 0: read N
                reader(buffer, &pop.N);
                unsigned nloci;
                reader(buffer, &nloci);
                // Read the mutations from the buffer
                io::read_mutations(buffer, pop.mutations);
                io::read_gametes(buffer, pop.gametes);
                unsigned ndips;
                reader(buffer, &ndips);
                pop.diploids.resize(
                    ndips, typename poptype::dipvector_t::value_type(nloci));
                io::deserialize_diploid<
                    typename poptype::dipvector_t::value_type::value_type>
                    dipreader;
                for (auto &dip : pop.diploids)
                    {
                        assert(dip.size() == nloci);
                        for (auto &genotype : dip)
                            {
                                dipreader(buffer, genotype);
                            }
                    }

                // update the mutation counts
                fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                                pop.mcounts);
                std::size_t temp;
                fwdpp::io::read_mutations(buffer, pop.fixations);
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
                            pop.mut_lookup.emplace(pop.mutations[i].pos,
                                                   static_cast<uint_t>(i));
                    }
            }
        }
    }
}
#endif
