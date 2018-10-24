#ifndef FWDPP_TS_DEFINITIONS_HPP
#define FWDPP_TS_DEFINITIONS_HPP

#include <cstdint>

/// \namespace fwdpp::ts Tree sequence \cite Kelleher2018-fu support

namespace fwdpp
{
    namespace ts
    {
        using TS_NODE_INT = std::int32_t;
        constexpr TS_NODE_INT TS_NULL_NODE = -1;
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
 * Note that this model has a rather big performance
 * bottleneck due to a suboptimal implementation of a 
 * "pick2" function.  We are working on it!
 */
#endif
