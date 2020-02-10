#ifndef FWDPP_TS_SIMPLIFICATION_ANCESTRY_LIST_HPP__
#define FWDPP_TS_SIMPLIFICATION_ANCESTRY_LIST_HPP__

#include <cstdint>
#include <algorithm>
#include <vector>
#include "segment.hpp"

namespace fwdpp
{
    namespace ts
    {
        class ancestry_list
        {
          private:
            void
            resize_and_fill(std::vector<std::int32_t>& v, std::size_t n)
            {
                v.resize(n);
                std::fill(begin(v), end(v), -1);
            }

          public:
            std::vector<segment> segments;
            std::vector<std::int32_t> first, next;

            ancestry_list() : segments(), first(), next() {}

            void
            init(std::size_t n)
            {
                segments.clear();
                resize_and_fill(first, n);
                resize_and_fill(next, n);
            }

            std::int32_t
            get_chain_tail(std::size_t i) const
            {
                if (i >= first.size())
                    {
                        throw std::runtime_error("index out of range");
                    }
                auto f = first[i];
                while (f != -1 && next[f] != -1)
                    {
                        f = next[f];
                    }
                return f;
            }

            void
            add_record(std::size_t i, double l, double r, TS_NODE_INT n)
            {
                if (i >= first.size())
                    {
                        throw std::runtime_error("index out of range");
                    }
                segments.emplace_back(l, r, n);
                if (first[i] == -1)
                    {
                        first[i] = segments.size() - 1;
                        if (segments.size() >= next.size())
                            {
                                next.push_back(-1);
                            }
                    }
                else
                    {
                        next.push_back(-1);
                        auto l = get_chain_tail(i);
                        next[l] = segments.size() - 1;
                    }
            }

            void
            nullify_chain(std::size_t i)
            {
                if (i >= first.size())
                    {
                        throw std::runtime_error("index out of range");
                    }
                first[i] = -1;
            }
        };

    } // namespace ts
} // namespace fwdpp

#endif
