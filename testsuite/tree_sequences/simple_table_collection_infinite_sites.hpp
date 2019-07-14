#ifndef FWDPP_TESTSUITE_SIMPLE_TABLE_COLLECTION_INFINITE_SITES_HPP
#define FWDPP_TESTSUITE_SIMPLE_TABLE_COLLECTION_INFINITE_SITES_HPP

#include <cmath>
#include "simple_table_collection.hpp"

class simple_table_collection_infinite_sites : public simple_table_collection
{
  private:
    const std::int8_t ancestral_state;
    const std::int8_t derived_state;

  public:
    simple_table_collection_infinite_sites()
        : simple_table_collection(), ancestral_state{ 0 }, derived_state{ 1 }
    {
        add_x_mutations_per_branch(1, 1);
    }

    void
    add_x_mutations_per_branch(int x, int exclude_root)
    {
        tables.site_table.clear();
        tables.mutation_table.clear();
        // 0.6 is arbitrary.  It just means we won't hit genome_length.
        double spacing_between_variants
            = 0.6 * tables.genome_length()
              / static_cast<double>(
                  x + (tables.node_table.size() - exclude_root));
        double next_variant_position = std::nexttoward(0., INFINITY);
        for (std::size_t i = 0; i < tables.node_table.size() - exclude_root;
             ++i)
            {
                for (int j = 0; j < x; ++j)
                    {
                        auto site = tables.emplace_back_site(
                            next_variant_position, ancestral_state);
                        tables.emplace_back_mutation(
                            static_cast<std::int32_t>(i), i + j, site,
                            derived_state, true);

                        next_variant_position += spacing_between_variants;
                    }
            }
    }
};

#endif

