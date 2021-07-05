#ifndef FWDPP_TS_MARGINAL_TREE_FUNCTIONS_SAMPLES_HPP__
#define FWDPP_TS_MARGINAL_TREE_FUNCTIONS_SAMPLES_HPP__

#include <stdexcept>
#include <fwdpp/util/named_type.hpp>
#include "../marginal_tree.hpp"
#include "../types/generate_null_id.hpp"

namespace fwdpp
{
    namespace ts
    {
        struct convert_sample_index_to_nodes_t
        {
        };

        /// Policy affecting behavior of samples_iterator
        using convert_sample_index_to_nodes
            = strong_types::named_type<bool, convert_sample_index_to_nodes_t>;

        template <typename SignedInteger> class samples_iterator
        /// \brief Faciliate traversal of the samples descending from a node
        /// \headerfile fwdpp/ts/marginal_tree_functions/samples.hpp
        {
          private:
            const marginal_tree<SignedInteger> &t;
            SignedInteger current_sample, right_sample;
            bool convert_to_nodes;

            inline const marginal_tree<SignedInteger> &
            init_marginal(const marginal_tree<SignedInteger> &m)
            {
                if (!m.advancing_sample_list())
                    {
                        throw std::invalid_argument(
                            "samples lists are not being updated");
                    }
                return m;
            }

            inline SignedInteger
            init_left_sample(SignedInteger u)
            {
                if (static_cast<std::size_t>(u) >= t.size())
                    {
                        throw std::invalid_argument("node index out of range");
                    }
                return t.left_sample[u];
            }

            inline SignedInteger
            init_right_sample(SignedInteger u)
            {
                if (static_cast<std::size_t>(u) >= t.size())
                    {
                        throw std::invalid_argument("node index out of range");
                    }
                return t.right_sample[u];
            }

          public:
            samples_iterator(const marginal_tree<SignedInteger> &m, SignedInteger u,
                             convert_sample_index_to_nodes convert)
                : t(init_marginal(m)), current_sample{init_left_sample(u)},
                  right_sample(init_right_sample(u)), convert_to_nodes(convert.get())
            /// \param m A marginal_tree
            /// \param u A node index
            {
            }

            inline SignedInteger
            operator()()
            /// Advance to the next sample
            ///
            /// Returns a null id when no more samples remain.
            {
                if (current_sample == types::generate_null_id<SignedInteger>())
                    {
                        //end of iteration
                        return current_sample;
                    }
                auto c = current_sample;
                if (c == right_sample)
                    {
                        // We are at the end of the samples list for this node,
                        // so we ensure that iteration will end
                        current_sample = types::generate_null_id<SignedInteger>();
                    }
                else
                    {
                        current_sample = t.next_sample[current_sample];
                    }
                return (convert_to_nodes) ? t.sample_table_index_to_node(c) : c;
            }

            template <typename F>
            inline bool
            operator()(const F &f)
            /// Apply a function to each sample
            /// \param f A function equivalent to void (*process_sample)(SignedInteger)
            ///
            /// Returns false to signify end of iteration.
            {
                auto s = this->operator()();
                bool rv = (s != types::generate_null_id<SignedInteger>());
                if (rv)
                    {
                        f(s);
                    }
                return rv;
            }
        };

        template <typename SignedInteger, typename F>
        inline void
        process_samples(const marginal_tree<SignedInteger> &m,
                        convert_sample_index_to_nodes convert, SignedInteger u,
                        const F &f)
        /// \brief Apply a function to the nodes descending from \a u
        /// \param m A fwdpp::ts::marginal_tree
        /// \param convert whether to iterate over sample nodes or sample list indexes
        /// \param u Index a fwdpp::ts::node
        /// \param f A function equivalent to void (*foo)(SignedInteger)
        {
            samples_iterator<SignedInteger> si(m, u, convert);
            while (si(f))
                {
                }
        }

        template <typename SignedInteger>
        inline std::vector<SignedInteger>
        get_samples(const marginal_tree<SignedInteger> &m, SignedInteger u)
        /// Return a vector of samples descending from node \a u in marginal_tree \a m.
        /// The return value contains node ids (as opposed to sample indexes)
        {
            std::vector<SignedInteger> rv;
            process_samples(m, convert_sample_index_to_nodes(true), u,
                            [&rv](SignedInteger x) { rv.push_back(x); });
            return rv;
        }

        template <typename SignedInteger>
        inline int
        num_samples(const marginal_tree<SignedInteger> &m, SignedInteger u)
        /// Return the number of samples descending from node \a u in marginal_tree \a m.
        {
            int nsamples = 0;
            process_samples(m, convert_sample_index_to_nodes(false), u,
                            [&nsamples](SignedInteger /*x*/) { ++nsamples; });
            return nsamples;
        }
    } // namespace ts
} // namespace fwdpp

#endif
