#ifndef FWDPP_IO_SERIALIZE_POPULATION_DETAIL_HPP__
#define FWDPP_IO_SERIALIZE_POPULATION_DETAIL_HPP__

#include <cstdint>
#include <stdexcept>
#include <fwdpp/io/scalar_serialization.hpp>
#include <fwdpp/io/mutation.hpp>
#include <fwdpp/io/haploid_genome.hpp>
#include <fwdpp/io/diploid.hpp>
#include <fwdpp/poptypes/tags.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

namespace fwdpp
{
    namespace io
    {
        namespace detail
        {
            template <typename streamtype, typename poptype>
            inline void
            serialize_population_details(streamtype &buffer, const poptype &pop,
                                         poptypes::DIPLOID_TAG)
            {
                io::scalar_writer writer;
                writer(buffer, &pop.N);
                io::write_mutations(buffer, pop.mutations);
                io::write_haploid_genomes(buffer, pop.haploid_genomes);
                io::write_diploids(buffer, pop.diploids);
                // Step 2: output fixations
                fwdpp::io::write_mutations(buffer, pop.fixations);
                if (!pop.fixations.empty())
                    {
                        // Step 3:the fixation times
                        writer(buffer, &pop.fixation_times[0], pop.fixations.size());
                    }
                // Write mcounts
                std::size_t num_mcounts = pop.mcounts.size();
                writer(buffer, &num_mcounts);
                if (num_mcounts)
                    {
                        writer(buffer, pop.mcounts.data(), pop.mcounts.size());
                    }
            }

            template <typename streamtype, typename poptype>
            inline void
            deserialize_population_details(poptype &pop, streamtype &buffer,
                                           poptypes::DIPLOID_TAG)
            {
                pop.clear();
                io::scalar_reader reader;
                // Step 0: read N
                reader(buffer, &pop.N);
                io::read_mutations(buffer, pop.mutations);
                io::read_haploid_genomes(buffer, pop.haploid_genomes);
                io::read_diploids(buffer, pop.diploids);

                // update the mutation counts
                // fwdpp_internal::process_haploid_genomes(pop.haploid_genomes,
                //                                        pop.mutations, pop.mcounts);
                fwdpp::io::read_mutations(buffer, pop.fixations);
                if (!pop.fixations.empty())
                    {
                        pop.fixation_times.resize(pop.fixations.size());
                        reader(buffer, &pop.fixation_times[0], pop.fixations.size());
                    }

                std::size_t num_mcounts{};
                reader(buffer, &num_mcounts);
                pop.mcounts.resize(num_mcounts);
                if (num_mcounts)
                    {
                        reader(buffer, pop.mcounts.data(), num_mcounts);
                    }

                // Finally, fill the lookup table:
                for (unsigned i = 0; i < pop.mcounts.size(); ++i)
                    {
                        if (pop.mcounts[i])
                            pop.mut_lookup.emplace(pop.mutations[i].pos,
                                                   static_cast<uint_t>(i));
                    }
            }
        } // namespace detail
    }     // namespace io
} // namespace fwdpp
#endif
