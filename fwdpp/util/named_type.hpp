#ifndef FWDPP_NAMED_TYPE_HPP_
#define FWDPP_NAMED_TYPE_HPP_

#include <type_traits>
#include <utility>

namespace fwdpp
{
    namespace strong_types
    {
        template <typename T, typename parameter_name> struct named_type
        /*! \brief Implementation of "strong types"
         *  
         *  Allows the construction of strong/unambiguous interfaces
         *  via thin aliases.  
         *
         *  At optimization levels like -O2 and higher, access 
         *  via get will be completely optimized out.
         *
         *  Example use includes fwdpp::trait and fwdpp::fitness
         *  as arguments to fwdpp::additive_diploid and
         *  fwdpp::multiplicative_diploid.
         *
         * \version 0.7.4 Added to fwdpp
         */
        {
          private:
            T value_;

          public:
            using value_type = T;
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
