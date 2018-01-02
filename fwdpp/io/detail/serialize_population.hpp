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
                                         sugar::SINGLEPOP_TAG)
            {
                io::scalar_writer writer;
                writer(buffer, &pop.N);
                io::write_mutations(pop.mutations, buffer);
                io::write_gametes(pop.gametes, buffer);
                io::write_diploids(pop.diploids, buffer);
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
            serialize_population_details(streamtype &buffer,
                                         const poptype &pop,
                                         sugar::MULTILOCPOP_TAG)
            {
                io::scalar_writer writer;
                writer(buffer, &pop.N);
                unsigned nloci = unsigned(pop.diploids[0].size());
                writer(buffer, &nloci);
                // write mutations
                io::write_mutations(pop.mutations, buffer);
                io::write_gametes(pop.gametes, buffer);
                unsigned ndips = unsigned(pop.diploids.size());
                writer(buffer, &ndips);
                io::serialize_diploid<
                    typename poptype::dipvector_t::value_type::value_type>
                    dipwriter;
                for (const auto &dip : pop.diploids)
                    {
                        for (const auto &genotype : dip)
                            {
                                dipwriter(genotype, buffer);
                            }
                    }
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
            serialize_population_details(streamtype &buffer,
                                         const poptype &pop,
                                         sugar::METAPOP_TAG)
            {
                io::scalar_writer writer;
                std::size_t npops = std::size_t(pop.Ns.size());
                writer(buffer, &npops);
                writer(buffer, &pop.Ns[0], npops);
                size_t i = unsigned(pop.diploids.size());
                writer(buffer, &i);
                io::write_mutations(pop.mutations, buffer);
                io::write_gametes(pop.gametes, buffer);
                for (const auto &deme : pop.diploids)
                    {
                        io::write_diploids(deme, buffer);
                    }
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
                io::read_mutations(pop.mutations, buffer);
                io::read_gametes(pop.gametes, buffer);
                io::read_diploids(pop.diploids, buffer);

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
                unsigned nloci;
                io::scalar_reader()(buffer, &nloci);
                // Read the mutations from the buffer
                io::read_mutations(pop.mutations, buffer);
                io::read_gametes(pop.gametes, buffer);
                unsigned ndips;
                io::scalar_reader()(buffer, &ndips);
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
                                dipreader(genotype, buffer);
                            }
                    }

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

                // Read diploids, mutaions, gametes
                std::size_t i;
                reader(buffer, & i);
                pop.diploids.resize(i);
                io::read_mutations(pop.mutations, buffer);
                io::read_gametes(pop.gametes, buffer);

                io::scalar_reader()(buffer, &i);
                for (auto &deme : pop.diploids)
                    {
                        io::read_diploids(deme, buffer);
                    }
                // update the mutation counts
                fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                                pop.mcounts);

                // read fixations
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
    }
}
#endif
