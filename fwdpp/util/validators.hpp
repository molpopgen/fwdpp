#ifndef FWDPP_UTIL_VALIDATORS_HPP
#define FWDPP_UTIL_VALIDATORS_HPP

#include <cmath>
#include <string>
#include <utility>
#include <stdexcept>
#include <type_traits>

namespace fwdpp
{
    namespace validators
    {
        template <typename E> struct unary_validator
        {
            template <typename T, typename F>
            inline T
            operator()(const T t, std::string message, const F f) const
            {
                if (!f(t))
                    {
                        throw E(std::move(message));
                    }
                return t;
            }
        };

        template <typename T>
        inline T
        casts_to_int(const T t, std::string message)
        {
            static_assert(std::is_arithmetic<T>::value,
                          "T must be an is_arithmetic type");
            static_assert(std::is_convertible<T, double>::value,
                          "T must be convertible to double");
            return unary_validator<std::invalid_argument>()(
                t, std::move(message), [](auto v) {
                    double p;
                    return std::modf(static_cast<double>(v), &p) == 0.;
                });
        }

        template <typename T>
        inline T
        non_negative(const T t, std::string message)
        {
            static_assert(std::is_arithmetic<T>::value,
                          "T must be an is_arithmetic type");
            return unary_validator<std::invalid_argument>()(
                t, std::move(message), [](auto v) { return v >= T(0); });
        }

        template <typename T>
        inline T
        is_positive(const T t, std::string message)
        {
            static_assert(std::is_arithmetic<T>::value,
                          "T must be an is_arithmetic type");
            return unary_validator<std::invalid_argument>()(
                t, std::move(message), [](auto v) { return v > T(0); });
        }

        template <typename T>
        inline T
        isfinite(const T t, std::string message)
        {
            static_assert(std::is_arithmetic<T>::value,
                          "T must be an is_arithmetic type");
            static_assert(std::is_convertible<T, double>::value,
                          "T must be convertible to double");
            return unary_validator<std::invalid_argument>()(
                t, std::move(message),
                [](auto v) { return std::isfinite(static_cast<double>(v)); });
        }

    }
}

#endif
