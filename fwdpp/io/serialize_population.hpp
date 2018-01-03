#ifndef FWDPP_IO_SERIALIZE_POPULATION_HPP__
#define FWDPP_IO_SERIALIZE_POPULATION_HPP__

#include "detail/serialize_population.hpp"

namespace fwdpp
{
    namespace io
    {
        template <typename streamtype, typename poptype>
        inline void
        serialize_population(streamtype &buffer, const poptype &pop)
        /// Write a population in binary format to a stream.
        ///
        /// \param buffer A model of std::ostream.
        /// \param pop A population.
        {
            detail::serialize_population_details(
                buffer, pop, typename poptype::popmodel_t());
        }

        template <typename streamtype, typename poptype>
        inline void
        deserialize_population(streamtype &buffer, poptype &pop)
        /// Read a binary representation of a population from a stream.
        ///
        /// \param pop A population to fill from \a buffer
        /// \param buffer A model of std::istream.
        {
            detail::deserialize_population_details(
                pop, buffer, typename poptype::popmodel_t());
        }
    }
}

#endif
