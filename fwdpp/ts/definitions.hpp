#ifndef FWDPP_TS_DEFINITIONS_HPP
#define FWDPP_TS_DEFINITIONS_HPP

#include <cstdint>

/// \namespace fwdpp::ts Tree sequence \cite Kelleher2018-fu support

namespace fwdpp
{
    namespace ts
    {
        /// Integer type for node indexes
        using TS_NODE_INT = std::int32_t;
        /// Index value of a NULL node
        constexpr TS_NODE_INT TS_NULL_NODE = -1;
        /// Convention for the ancestral state of a site
        constexpr std::int8_t default_ancestral_state = 0;
        /// Convention for the derived state of a site
        constexpr std::int8_t default_derived_state = 1;
    } // namespace ts
} // namespace fwdpp

/*! \example wfts.cc 
*
*  Example of Wright-Fisher simulation with tree sequences
*/

/*! \example spatialts.cc
 *
 * Example of continous space and tree sequences.
 *
 * This example uses simple methods to find mates
 * within a euclidean distance of an individual, and
 * then choose a mate proportional to fitnesses within
 * that circle.
 *
 * Note that this example has a rather big performance
 * bottleneck due to a suboptimal implementation of a 
 * "pick2" function.  We are working on it!
 */
#endif
