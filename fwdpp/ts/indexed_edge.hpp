#ifndef FWDPP_TS_INDEXED_EDGE_HPP
#define FWDPP_TS_INDEXED_EDGE_HPP

#include <cstdint>
#include <vector>

namespace fwdpp
{
    namespace ts
    {
        struct indexed_edge
        /// Holds edge data keyed on either right or left position.
        /// Used to define the index vectors described on page 13
        /// of \cite Kelleher2016-cb
        ///  \version 0.7.0 Added to fwdpp
        {
            double pos, time;
            std::int32_t parent, child;
            explicit indexed_edge(double pos_, double t, std::int32_t p,
                                  std::int32_t c)
                : pos{ pos_ }, time{ t }, parent{ p }, child{ c }
            {
            }
            inline bool
            operator<(const indexed_edge& rhs) const
            {
                if (pos == rhs.pos)
                    {
                        return time < rhs.time;
                    }
                return pos < rhs.pos;
            }
        };

		/// An index for an edge table.  See \cite Kelleher2016-cb, 
		/// page 13
        ///  \version 0.7.0 Added to fwdpp
        using indexed_edge_container = std::vector<indexed_edge>;
    } // namespace ts
} // namespace fwdpp

#endif
