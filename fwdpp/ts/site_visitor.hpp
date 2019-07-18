#ifndef FWDPP_TS_SITE_VISITOR_HPP
#define FWDPP_TS_SITE_VISITOR_HPP

#include <algorithm>
#include <utility>
#include "tree_visitor.hpp"

namespace fwdpp
{
    namespace ts
    {
        class site_visitor
        /// \brief Facilitae iteration over ts::site objects one at a time
        /// For example use, see implementation of ts::generate_data_matrix.
        /// \version 0.8.0 Added to fwdpp
        {
          private:
            const table_collection& tables_;
            tree_visitor tv;
            site_vector::const_iterator current_site;
            mutation_key_vector::const_iterator current_mutation;
            std::pair<mutation_key_vector::const_iterator,
                      mutation_key_vector::const_iterator>
                mutations_at_current_site;

            void
            set_mutations_at_current_site(
                site_vector::const_iterator itr,
                mutation_key_vector::const_iterator mitr)
            {
                if (itr == std::end(tables_.site_table))
                    {
                        mutations_at_current_site
                            = { std::end(tables_.mutation_table),
                                std::end(tables_.mutation_table) };
                    }
                auto m = std::find_if(
                    mitr, std::end(tables_.mutation_table),
                    [this, itr](const mutation_record& mr) {
                        return tables_.site_table[mr.site].position
                               != itr->position;
                    });
                mutations_at_current_site = { mitr, m };
            }

            template <typename SAMPLES>
            tree_visitor
            init_tree_visitor(const SAMPLES& samples)
            {
                tree_visitor tv(tables_, samples, update_samples_list(true));
                auto t = tv();
                if (!t)
                    {
                        throw tables_error("no tree in table collection");
                    }
                return tv;
            }

          public:
            template <typename SAMPLES>
            site_visitor(const table_collection& tables,
                         const SAMPLES& samples)
                : tables_(tables), tv(init_tree_visitor(samples)),
                  current_site(begin(tables.site_table)),
                  current_mutation(begin(tables.mutation_table)),
                  mutations_at_current_site(std::end(tables_.mutation_table),
                                            std::end(tables_.mutation_table))
            {
            }

            site_vector::const_iterator
            operator()()
            {
                if (current_site == std::end(tables_.site_table))
                    {
                        return current_site;
                    }
                // Need to advance trees here if needed
                while (current_site < std::end(tables_.site_table)
                       && (current_site->position < tv.tree().left
                           || current_site->position >= tv.tree().right))
                    {
                        auto t = tv();
                        if (!t && current_site < std::end(tables_.site_table))
                            {
                                throw std::runtime_error(
                                    "tree sequence interation error");
                            }
                    }
                while (current_mutation < std::end(tables_.mutation_table)
                       && tables_.site_table[current_mutation->site].position
                              < current_site->position)
                    {
                        ++current_mutation;
                    }
                if (current_mutation < std::end(tables_.mutation_table)
                    && tables_.site_table[current_mutation->site].position
                           != current_site->position)
                    {
                        throw tables_error(
                            "site and mutation tables are invalid");
                    }
                auto temp = current_site;
                set_mutations_at_current_site(current_site, current_mutation);
                ++current_site;
                return temp;
            }

            site_vector::const_iterator
            end() const
            {
                return std::end(tables_.site_table);
            }

            std::pair<mutation_key_vector::const_iterator,
                      mutation_key_vector::const_iterator>
            get_mutations() const
            {
                return mutations_at_current_site;
            }

            const marginal_tree&
            current_tree() const
            {
                return tv.tree();
            }
        };

        inline site_vector::const_iterator
        end(site_visitor& sv)
        {
            return sv.end();
        }
    } // namespace ts
} // namespace fwdpp

#endif

