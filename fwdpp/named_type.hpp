#ifndef FWDPP_NAMED_TYPE_HPP_
#define FWDPP_NAMED_TYPE_HPP_

#include <type_traits>
#include <utility>

namespace fwdpp
{
    namespace strong_types
    {
        template <typename T, typename parameter_name> struct named_type
        {
          private:
            T value_;

          public:
            explicit named_type(const T& value) : value_(value) {}
            template <typename T_>
            explicit named_type(
                T_&& value,
                typename std::enable_if<!std::is_reference<T_>::value>::value)
                : value_(std::move(value))
            {
            }

            T&
            get()
            {
                return value_;
            }

            const T&
            get() const
            {
                return value_;
            }
        };
    } // namespace strong_types
} // namespace fwdpp

#endif
