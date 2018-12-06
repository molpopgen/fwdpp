#ifndef FWDPP_WRAPPED_RANGE_HPP
#define FWDPP_WRAPPED_RANGE_HPP

#include <utility>

namespace fwdpp
{
    template <typename T> struct wrapped_range
    {
        T begin;
        T end;
        template <typename T_>
        wrapped_range(T_ &&b, T_ &&e)
            : begin{ std::forward<T_>(b) }, end{ std::forward<T_>(e) }
        {
        }
    };

    template <typename T>
    inline wrapped_range<T>
    make_wrapped_range(T b, T e)
    {
        return wrapped_range<T>(b, e);
    }

    template <typename T>
    inline T
    begin(wrapped_range<T> &wr)
    {
        return wr.begin;
    }

    template <typename T>
    inline auto
    begin(const wrapped_range<T> &wr) -> decltype(wr.begin)
    {
        return wr.begin;
    }

    template <typename T>
    inline T
    end(wrapped_range<T> &wr)
    {
        return wr.end;
    }

    template <typename T>
    inline auto
    end(const wrapped_range<T> &wr) -> decltype(wr.end)
    {
        return wr.end;
    }
} // namespace fwdpp

#endif
