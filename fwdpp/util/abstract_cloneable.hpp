#ifndef FWDPP_UTIL_ABSTRACT_CLONEABLE_HPP__
#define FWDPP_UTIL_ABSTRACT_CLONEABLE_HPP__

#include <memory>

namespace fwdpp
{
    namespace util
    {
        template <typename T> struct abstract_cloneable
        /// Generate cloneable abstract base classes via CRTP
        {
			/// Return a clone (new copy) of object
            virtual std::unique_ptr<T> clone() const = 0;
            abstract_cloneable() = default;
            virtual ~abstract_cloneable() = default;
            abstract_cloneable(const abstract_cloneable &) = delete;
            abstract_cloneable(abstract_cloneable &&) = delete;
            abstract_cloneable &operator=(const abstract_cloneable &) = delete;
            abstract_cloneable &operator=(abstract_cloneable &&) = delete;
        };
    } // namespace util
} // namespace fwdpp

#endif
