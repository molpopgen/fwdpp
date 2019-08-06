#ifndef FWDPP_TS_DETAIL_SPLIT_BREAKPOINTS_HPP
#define FWDPP_TS_DETAIL_SPLIT_BREAKPOINTS_HPP

#include <vector>
#include <limits>
#include <algorithm>
#include <cassert>
#include <tuple>
#include "../definitions.hpp"
#include "../edge.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace detail
        {
            template <typename F>
            void
            split_breakpoints_add_edges(
                const std::vector<double>& breakpoints,
                const std::tuple<TS_NODE_INT, TS_NODE_INT>& parents,
                const TS_NODE_INT next_index, const double genome_length,
                const F& f)
            {
                if (breakpoints.front() != 0.0)
                    {
                        f(edge{ 0.0, breakpoints.front(), std::get<0>(parents),
                                next_index });
                        //this->push_back_edge(0., breakpoints.front(),
                        //                     std::get<0>(parents), next_index);
                    }
                // TODO: replace with exception via a debug mode
                assert(std::count(begin(breakpoints), end(breakpoints),
                                  std::numeric_limits<double>::max())
                       == 1);
                assert(breakpoints.back()
                       == std::numeric_limits<double>::max());
                for (unsigned j = 1; j < breakpoints.size(); ++j)
                    {
                        double a = breakpoints[j - 1];
                        double b = (j < breakpoints.size() - 1)
                                       ? breakpoints[j]
                                       : genome_length;
                        if (b <= a)
                            {
                                throw std::invalid_argument(
                                    "right must be > left");
                            }
                        if (j % 2 == 0.)
                            {
                                f(edge{ a, b, std::get<0>(parents),
                                        next_index });
                            }
                        else
                            {
                                f(edge{ a, b, std::get<1>(parents),
                                        next_index });
                            }
                    }
            }

            template <typename F>
            void
            split_breakpoints(
                const std::vector<double>& breakpoints,
                const std::tuple<TS_NODE_INT, TS_NODE_INT>& parents,
                const TS_NODE_INT next_index, const double genome_length,
                const F& f)
            {
                if (breakpoints.empty())
                    {
                        f(edge{ 0., genome_length, std::get<0>(parents),
                                next_index });
                        return;
                    }
                auto itr = std::adjacent_find(std::begin(breakpoints),
                                              std::end(breakpoints));
                if (itr == std::end(breakpoints))
                    {
                        split_breakpoints_add_edges(breakpoints, parents,
                                                    next_index, genome_length,
                                                    f);
                    }
                else
                    {
                        // Here, we need to reduce the input
                        // breakpoints to only those seen
                        // an odd number of times.
                        // Even numbers of the same breakpoint
                        // are "double x-overs" and thus
                        // cannot affect the genealogy.
                        std::vector<double> odd_breakpoints;
                        auto start = breakpoints.begin();
                        while (itr < breakpoints.end())
                            {
                                auto not_equal
                                    = std::find_if(itr, breakpoints.end(),
                                                   [itr](const double d) {
                                                       return d != *itr;
                                                   });
                                int even = (std::distance(itr, not_equal) % 2
                                            == 0.0);
                                odd_breakpoints.insert(odd_breakpoints.end(),
                                                       start, itr + 1 - even);
                                start = not_equal;
                                itr = std::adjacent_find(
                                    start, std::end(breakpoints));
                            }
                        odd_breakpoints.insert(odd_breakpoints.end(), start,
                                               breakpoints.end());
                        split_breakpoints_add_edges(odd_breakpoints, parents,
                                                    next_index, genome_length,
                                                    f);
                    }
            }
        } // namespace detail
    }     // namespace ts
} // namespace fwdpp

#endif
