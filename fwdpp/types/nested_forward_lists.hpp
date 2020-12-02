#ifndef FWDPP_TYPES_NESTED_FORWARD_LISTS_HPP
#define FWDPP_TYPES_NESTED_FORWARD_LISTS_HPP

#include <vector>
#include <limits>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <string>
#include <utility>

namespace fwdpp
{
    class nested_forward_lists_overflow : public std::exception
    {
      private:
        std::string message_;

      public:
        explicit nested_forward_lists_overflow(std::string message)
            : message_(std::move(message))
        {
        }
        virtual const char*
        what() const noexcept
        {
            return message_.c_str();
        }
    };

    template <typename T, typename Index, Index NullValue> class nested_forward_lists
    /// Container of multiple owning forward lists with corresponding head/tail index vector
    /// describing where individual lists start/stop.
    /// Addition of a list at a new index creates null entries for intervening head/tail
    /// indexes as needed.
    /// Intended use is to track linked lists whose entry points (head) can be specified
    /// by integers.
    /// Const forward/backward iterator access is provided to the head vector.
    ///
    /// \version 0.9.0 Added to library
    {
      private:
        static_assert(std::is_integral<Index>::value, "Index must be an integer type");

        template <typename... Args>
        void
        insert_new_record(std::size_t idx, Args&&... args)
        {
            data.emplace_back(std::forward<Args>(args)...);
            _head[idx] = static_cast<Index>(data.size() - 1);
            _tail[idx] = _head[idx];
            _next.emplace_back(null);
        }

        void
        throw_if_null(Index i) const
        {
            if (i == null)
                {
                    throw std::invalid_argument("index is null");
                }
        }

        void
        validate_index(Index i, std::size_t size) const
        {
            if (static_cast<std::size_t>(i) >= size)
                {
                    throw std::out_of_range("index out of range");
                }
        }

        template <typename Container>
        void
        swap_with_empty(Container& c)
        {
            Container temp;
            c.swap(temp);
        }

        std::vector<T> data;
        std::vector<Index> _head, _tail, _next;

      public:
        static constexpr Index null = NullValue;

        using const_iterator = typename std::vector<Index>::const_iterator;
        using const_reverse_iterator =
            typename std::vector<Index>::const_reverse_iterator;

        nested_forward_lists() : data{}, _head{}, _tail{}, _next{}
        {
        }

        Index
        tail(Index at) const
        {
            throw_if_null(at);
            validate_index(at, _tail.size());
            return _tail[static_cast<std::size_t>(at)];
        }

        Index
        head(Index at) const
        {
            throw_if_null(at);
            validate_index(at, _head.size());
            return _head[static_cast<std::size_t>(at)];
        }

        Index
        next(Index at) const
        {
            throw_if_null(at);
            validate_index(at, data.size());
            return _next[static_cast<std::size_t>(at)];
        }

        template <typename... Args>
        void
        extend(Index at, Args&&... args)
        {
            throw_if_null(at);
            if (data.size() >= std::numeric_limits<Index>::max() - 1)
                {
                    throw nested_forward_lists_overflow(
                        "buffer has overflowed Index maximum");
                }
            auto idx = static_cast<std::size_t>(at);
            if (idx >= _head.size())
                {
                    _head.resize(idx + 1, null);
                    _tail.resize(idx + 1, null);
                }
            if (_head[idx] == null)
                {
                    insert_new_record(idx, std::forward<Args>(args)...);
                    return;
                }
            auto t = _tail[idx];
            if (t == null)
                {
                    throw std::runtime_error("unexpected null tail value");
                }
            data.emplace_back(std::forward<Args>(args)...);
            _tail[idx] = static_cast<Index>(data.size() - 1);
            _next[static_cast<std::size_t>(t)] = static_cast<Index>(data.size() - 1);
            _next.emplace_back(null);
        }

        template <typename Container>
        void
        extend_from_container(Index at, const Container& c)
        {
            for (auto& i : c)
                {
                    extend(at, i);
                }
        }

        T&
        fetch(Index at)
        {
            throw_if_null(at);
            validate_index(at, data.size());
            return data[static_cast<std::size_t>(at)];
        }

        const T&
        fetch(Index at) const
        {
            throw_if_null(at);
            validate_index(at, data.size());
            return data[static_cast<std::size_t>(at)];
        }

        void
        nullify_list(Index at)
        // In future, we can recycle indexes in the
        // data vector from this list for re-use.
        {
            throw_if_null(at);
            validate_index(at, _head.size());
            auto idx = static_cast<std::size_t>(at);
            _head[idx] = _tail[idx] = null;
        }

        void
        reset(std::size_t newsize)
        {
            clear();
            _head.resize(newsize);
            _tail.resize(newsize);
            std::fill(std::begin(_head), std::end(_head), null);
            std::fill(std::begin(_tail), std::end(_tail), null);
        }

        void
        clear()
        {
            data.clear();
            _head.clear();
            _tail.clear();
            _next.clear();
        }

        void
        release_memory()
        {
            swap_with_empty(data);
            swap_with_empty(_head);
            swap_with_empty(_tail);
            swap_with_empty(_next);
        }

        const_iterator
        begin() const
        {
            return _head.begin();
        }

        const_iterator
        end() const
        {
            return _head.end();
        }

        const_reverse_iterator
        rbegin() const
        {
            return _head.rbegin();
        }

        const_reverse_iterator
        rend() const
        {
            return _head.rend();
        }

        Index
        convert_to_head_index(const_iterator i) const
        {
            return std::distance(_head.cbegin(), i);
        }

        Index
        convert_to_head_index(const_reverse_iterator i) const
        {
            return convert_to_head_index(i.base() - 1);
        }
    };

    template <typename T, typename Index, Index NullValue>
    constexpr Index nested_forward_lists<T, Index, NullValue>::null;

    //template <typename nested_forward_lists_t>
    //inline typename nested_forward_lists_t::const_iterator
    //begin(const nested_forward_lists_t& n)
    //{
    //    return n.begin();
    //}

    //template <typename nested_forward_lists_t>
    //inline typename nested_forward_lists_t::const_iterator
    //end(const nested_forward_lists_t& n)
    //{
    //    return n.end();
    //}

    template <typename nested_forward_lists_t>
    inline typename nested_forward_lists_t::const_reverse_iterator
    rbegin(const nested_forward_lists_t& n)
    {
        return n.rbegin();
    }

    template <typename nested_forward_lists_t>
    inline typename nested_forward_lists_t::const_reverse_iterator
    rend(const nested_forward_lists_t& n)
    {
        return n.rend();
    }
}
#endif

