#ifndef FWDPP_TS_SITE_VISITOR_HPP
#define FWDPP_TS_SITE_VISITOR_HPP

#include <algorithm>
#include <utility>
#include "tree_visitor.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename TableCollectionType> class site_visitor
        /// \brief Facilitae iteration over ts::site objects one at a time
        /// For example use, see implementation of ts::generate_data_matrix.
        /// \version 0.8.0 Added to fwdpp
        /// \version 0.9.0 Made a template class
        {
          private:
            const TableCollectionType& tables_;
            tree_visitor<TableCollectionType> tv;
            typename TableCollectionType::site_table::const_iterator current_site;
            typename TableCollectionType::mutation_table::const_iterator
                current_mutation;
            std::pair<typename TableCollectionType::mutation_table::const_iterator,
                      typename TableCollectionType::mutation_table::const_iterator>
                mutations_at_current_site;

            void
            set_mutations_at_current_site(
                typename TableCollectionType::site_table::const_iterator itr,
                typename TableCollectionType::mutation_table::const_iterator mitr)
            {
                if (itr == std::end(tables_.sites))
                    {
                        mutations_at_current_site
                            = {std::end(tables_.mutations), std::end(tables_.mutations)};
                    }
                auto m = std::find_if(
                    mitr, std::end(tables_.mutations),
                    [this,
                     itr](typename TableCollectionType::mutation_table::const_reference
                              mr) {
                        return tables_.sites[mr.site].position != itr->position;
                    });
                mutations_at_current_site = {mitr, m};
            }

            template <typename SAMPLES>
            tree_visitor<TableCollectionType>
            init_tree_visitor(const SAMPLES& samples)
            {
                tree_visitor<TableCollectionType> tv(tables_, samples, update_samples_list(true));
                auto t = tv();
                if (!t)
                    {
                        throw tables_error("no tree in table collection");
                    }
                return tv;
            }

          public:
            template <typename SAMPLES>
            site_visitor(const TableCollectionType& tables, const SAMPLES& samples)
                : tables_(tables), tv(init_tree_visitor(samples)),
                  current_site(begin(tables.sites)),
                  current_mutation(begin(tables.mutations)),
                  mutations_at_current_site(std::end(tables_.mutations),
                                            std::end(tables_.mutations))
            {
            }

            typename TableCollectionType::site_table::const_iterator
            operator()()
            {
                if (current_site == std::end(tables_.sites))
                    {
                        return current_site;
                    }
                // Need to advance trees here if needed
                while (current_site < std::end(tables_.sites)
                       && (current_site->position < tv.tree().left
                           || current_site->position >= tv.tree().right))
                    {
                        auto t = tv();
                        if (!t && current_site < std::end(tables_.sites))
                            {
                                throw std::runtime_error(
                                    "tree sequence interation error");
                            }
                    }
                while (current_mutation < std::end(tables_.mutations)
                       && tables_.sites[current_mutation->site].position
                              < current_site->position)
                    {
                        ++current_mutation;
                    }
                if (current_mutation < std::end(tables_.mutations)
                    && tables_.sites[current_mutation->site].position
                           != current_site->position)
                    {
                        throw tables_error("site and mutation tables are invalid");
                    }
                auto temp = current_site;
                set_mutations_at_current_site(current_site, current_mutation);
                ++current_site;
                return temp;
            }

            typename TableCollectionType::site_table::const_iterator
            end() const
            {
                return std::end(tables_.sites);
            }

            std::pair<typename TableCollectionType::mutation_table::const_iterator,
                      typename TableCollectionType::mutation_table::const_iterator>
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

        template <typename TableCollectionType>
        inline typename TableCollectionType::site_table_t::const_iterator
        end(site_visitor<TableCollectionType>& sv)
        {
            return sv.end();
        }
    } // namespace ts
} // namespace fwdpp

#endif

