#ifndef __FWDPP_SUGAR_SINGLEPOP_SERIALIZATON_HPP__
#define __FWDPP_SUGAR_SINGLEPOP_SERIALIZATON_HPP__

#include <ios>
#include <iosfwd>
#include <sstream>
#include <algorithm>
#include <type_traits>
#include <utility>
#include <functional>
#include <fwdpp/IO.hpp>
#include <fwdpp/sugar/poptypes/tags.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/generalmut.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>

namespace KTfwd
{
    /*!
      \brief Facilitates serialization of mutation
      types supported by the fwdpp sugar library
      \ingroup sugar
    */
    struct mutation_writer
    {
        /*!
          \brief overload for KTfwd::popgenmut and ostreams
         */
        using result_type = void;
        template <typename streamtype, typename mutation_t>
        inline
            typename std::enable_if<std::is_same<mutation_t, popgenmut>::value,
                                    result_type>::type
            operator()(const mutation_t &m, streamtype &buffer) const
        {
            fwdpp_internal::scalar_writer writer;
            writer(buffer, &m.g);
            writer(buffer, &m.pos);
            writer(buffer, &m.s);
            writer(buffer, &m.h);
        }

        template <typename streamtype, typename mutation_t>
        inline
            typename std::enable_if<std::is_same<mutation_t, mutation>::value,
                                    result_type>::type
            operator()(const mutation_t &m, streamtype &buffer) const
        {
            fwdpp_internal::scalar_writer writer;
            writer(buffer, &m.pos);
            writer(buffer, &m.s);
            writer(buffer, &m.h);
        }
        //! \brief overload for KTfwd::generalmut and ostream
        template <typename streamtype, typename mutation_t,
                  std::size_t N
                  = std::tuple_size<typename mutation_t::array_t>::value>
        inline typename std::enable_if<std::is_same<mutation_t,
                                                    generalmut<N>>::value,
                                       result_type>::type
        operator()(const mutation_t &t, streamtype &buffer) const
        {
            fwdpp_internal::scalar_writer writer;
            writer(buffer, &t.g);
            writer(buffer, &t.pos);
            // Write mutation types
            writer(buffer, &t.s[0], N);
            writer(buffer, &t.h[0], N);
        }

        //! \brief overload for KTfwd::generalmut_vec and ostream
        template <typename streamtype, typename mutation_t>
        inline typename std::enable_if<std::is_same<mutation_t,
                                                    generalmut_vec>::value,
                                       result_type>::type
        operator()(const mutation_t &t, streamtype &buffer) const
        {
            fwdpp_internal::scalar_writer writer;
            writer(buffer, &t.g);
            writer(buffer, &t.pos);
            // Write mutation types
            using array_t_size_t = typename generalmut_vec::array_t::size_type;
            array_t_size_t ns = t.s.size(), nh = t.h.size();
            writer(buffer, &ns);
            writer(buffer, &nh);
            writer(buffer, t.s.data(), ns);
            writer(buffer, t.h.data(), nh);
        }
    };
    /*!
      \brief Facilitates serialization of mutation
      types supported by the fwdpp sugar library.

      The template parameter must be derived from
      KTfwd::mutation_base

      \ingroup sugar
    */
    template <typename mutation_t> struct mutation_reader
    {
        //! The return value of operator()
        using result_type = mutation_t;
        /*!
          \brief overload for KTfwd::popgenmut and istreams
         */
        template <typename streamtype, typename U = mutation_t>
        inline typename std::enable_if<std::is_same<U, popgenmut>::value,
                                       result_type>::type
        operator()(streamtype &in) const
        {
            uint_t g;
            double pos, s, h;
            fwdpp_internal::scalar_reader reader;
            reader(in, &g);
            reader(in, &pos);
            reader(in, &s);
            reader(in, &h);
            return result_type(pos, s, h, g);
        }
        /*!
          \brief overload for KTfwd::mutation and istreams
         */
        template <typename streamtype, typename U = mutation_t>
        inline typename std::enable_if<std::is_same<U, mutation>::value,
                                       result_type>::type
        operator()(streamtype &in) const
        {
            double pos, s, h;
            fwdpp_internal::scalar_reader reader;
            reader(in, &pos);
            reader(in, &s);
            reader(in, &h);
            return result_type(pos, s, h);
        }
        //! \brief overalod for KTfwd::generalmut and std::istream
        template <typename streamtype, typename U = mutation_t>
        inline typename std::
            enable_if<std::is_same<U, generalmut<std::tuple_size<
                                          typename U::array_t>::value>>::value,
                      result_type>::type
            operator()(streamtype &in) const
        {
            uint_t g;
            double pos;
            using value_t = typename U::array_t::value_type;
            fwdpp_internal::scalar_reader reader;
            std::array<value_t, std::tuple_size<typename U::array_t>::value> s,
                h;
            reader(in, &g);
            reader(in, &pos);
            // Write mutation types
            reader(in, &s[0], std::tuple_size<typename U::array_t>::value);
            reader(in, &h[0], std::tuple_size<typename U::array_t>::value);
            return generalmut<std::tuple_size<typename U::array_t>::value>(
                s, h, pos, g);
        }

        //! \brief overalod for KTfwd::generalmut_vec and std::istream
        template <typename streamtype, typename U = mutation_t>
        inline typename std::enable_if<std::is_same<U, generalmut_vec>::value,
                                       result_type>::type
        operator()(streamtype &in) const
        {
            uint_t g;
            double pos;
            fwdpp_internal::scalar_reader reader;
            reader(in, &g);
            reader(in, &pos);
            typename U::array_t::size_type ns, nh;
            reader(in, &ns);
            reader(in, &nh);
            typename U::array_t s(ns), h(nh);
            // Write mutation types
            if (ns)
                {
                    reader(in, s.data(), ns);
                }
            if (nh)
                {
                    reader(in, h.data(), nh);
                }
            return U(std::move(s), std::move(h), pos, g);
        }
    };

    /*!
      \brief Serialize populations.
      \ingroup sugar
     */
    struct serialize
    {
        using result_type = void;
        fwdpp_internal::scalar_writer writer;
        //! Default constructor
        serialize() : writer(fwdpp_internal::scalar_writer()) {}

        //! Move constructor.  Req'd for this to be a member type of another
        //! class
        serialize(serialize &&) = default;
        //    serialize( serialize && __s) : buffer(__buffer.str()) {
        //    }

        /*!
          \brief Overload for single population simulations
         */
        template <typename streamtype, typename sugarpop_t, typename writer_t,
                  typename diploid_writer_t = diploidIOplaceholder>
        inline typename std::
            enable_if<std::is_same<typename sugarpop_t::popmodel_t,
                                   sugar::SINGLEPOP_TAG>::value,
                      result_type>::type
            operator()(streamtype &buffer, const sugarpop_t &pop,
                       const writer_t &wt,
                       const diploid_writer_t &dw = diploid_writer_t()) const
        {
            writer(buffer, &pop.N);
            write_binary_pop(pop.gametes, pop.mutations, pop.diploids, wt,
                             buffer, dw);
            // Step 2: output fixations
            uint_t temp = uint_t(pop.fixations.size());
            writer(buffer, &temp);
            if (temp)
                {
                    std::for_each(pop.fixations.begin(), pop.fixations.end(),
                                  std::bind(std::cref(wt),
                                            std::placeholders::_1,
                                            std::ref(buffer)));
                    // Step 3:the fixation times
                    writer(buffer, &pop.fixation_times[0], temp);
                }
        }

        template <typename streamtype, typename sugarpop_t, typename writer_t,
                  typename diploid_writer_t = diploidIOplaceholder>
        inline typename std::
            enable_if<std::is_same<typename sugarpop_t::popmodel_t,
                                   sugar::MULTILOCPOP_TAG>::value,
                      result_type>::type
            operator()(streamtype &buffer, const sugarpop_t &pop,
                       const writer_t &wt,
                       const diploid_writer_t &dw = diploid_writer_t()) const
        {
            writer(buffer, &pop.N);
            write_binary_pop_mloc(pop.gametes, pop.mutations, pop.diploids, wt,
                                  buffer, dw);
            // Step 2: output fixations
            uint_t temp = uint_t(pop.fixations.size());
            writer(buffer, &temp);
            if (temp)
                {
                    std::for_each(pop.fixations.begin(), pop.fixations.end(),
                                  std::bind(std::cref(wt),
                                            std::placeholders::_1,
                                            std::ref(buffer)));
                    // Step 3:the fixation times
                    writer(buffer, &pop.fixation_times[0], temp);
                }
            temp = uint_t(pop.locus_boundaries.size());
            writer(buffer, &temp);
            for(uint_t i=0;i<temp;++i)
            {
                writer(buffer,&pop.locus_boundaries[i].first);
                writer(buffer,&pop.locus_boundaries[i].second);
            }
        }

        /*!
          \brief Overload for metapopulation simulations
         */
        template <typename streamtype, typename sugarpop_t, typename writer_t,
                  typename diploid_writer_t = diploidIOplaceholder>
        inline typename std::
            enable_if<std::is_same<typename sugarpop_t::popmodel_t,
                                   sugar::METAPOP_TAG>::value,
                      result_type>::type
            operator()(streamtype &buffer, const sugarpop_t &pop,
                       const writer_t &wt,
                       const diploid_writer_t &dw = diploid_writer_t()) const
        {
            uint_t npops = uint_t(pop.Ns.size());
            writer(buffer, &npops);
            writer(buffer, &pop.Ns[0], npops);
            write_binary_metapop(pop.gametes, pop.mutations, pop.diploids, wt,
                                 buffer, dw);
            // Step 2: output fixations
            uint_t temp = uint_t(pop.fixations.size());
            writer(buffer, &temp);
            if (temp)
                {
                    std::for_each(pop.fixations.begin(), pop.fixations.end(),
                                  std::bind(std::cref(wt),
                                            std::placeholders::_1,
                                            std::ref(buffer)));
                    // Step 3:the fixation times
                    writer(buffer, &pop.fixation_times[0], temp);
                }
        }
    };

    /*!
      \brief Deserialize population objects
      \ingroup sugar
     */
    struct deserialize
    {
        //! The return type for operator()
        using result_type = void;
        /*!
          \brief Overload for single population simulations
         */
        template <typename streamtype, typename sugarpop_t, typename reader_t,
                  typename diploid_reader_t = diploidIOplaceholder>
        inline typename std::
            enable_if<std::is_same<typename sugarpop_t::popmodel_t,
                                   sugar::SINGLEPOP_TAG>::value,
                      result_type>::type
            operator()(sugarpop_t &pop, streamtype &buffer, const reader_t &rt,
                       const diploid_reader_t &dr = diploid_reader_t()) const
        {
            pop.clear();
            fwdpp_internal::scalar_reader reader;
            // Step 0: read N
            reader(buffer, &pop.N);
            KTfwd::read_binary_pop(pop.gametes, pop.mutations, pop.diploids,
                                   rt, buffer, dr);

            // update the mutation counts
            fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                            pop.mcounts);
            uint_t temp;
            reader(buffer, &temp);
            for (uint_t m = 0; m < temp; ++m)
                {
                    typename reader_t::result_type mm = rt(buffer);
                    pop.fixations.emplace_back(std::move(mm));
                }
            pop.fixation_times.resize(temp);
            if (temp)
                {
                    reader(buffer, &pop.fixation_times[0], temp);
                }

            // Finally, fill the lookup table:
            for (unsigned i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (pop.mcounts[i])
                        pop.mut_lookup.insert(pop.mutations[i].pos);
                }
        }

        template <typename streamtype, typename sugarpop_t, typename reader_t,
                  typename diploid_reader_t = diploidIOplaceholder>
        inline typename std::
            enable_if<std::is_same<typename sugarpop_t::popmodel_t,
                                   sugar::MULTILOCPOP_TAG>::value,
                      result_type>::type
            operator()(sugarpop_t &pop, streamtype &buffer, const reader_t &rt,
                       const diploid_reader_t &dr = diploid_reader_t()) const
        {
            pop.clear();
            fwdpp_internal::scalar_reader reader;
            // Step 0: read N
            reader(buffer, &pop.N);
            KTfwd::read_binary_pop_mloc(pop.gametes, pop.mutations,
                                        pop.diploids, rt, buffer, dr);
            // update the mutation counts
            fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                            pop.mcounts);
            uint_t temp;
            reader(buffer, &temp);
            for (uint_t m = 0; m < temp; ++m)
                {
                    typename reader_t::result_type mm = rt(buffer);
                    pop.fixations.emplace_back(std::move(mm));
                }
            pop.fixation_times.resize(temp);
            if (temp)
                {
                    reader(buffer, &pop.fixation_times[0], temp);
                }

            reader(buffer,&temp);
            if(temp)
            {
                double x[2];
                for(uint_t i=0;i<temp;++i)
                {
                    reader(buffer,&x[0],2);
                    pop.locus_boundaries.emplace_back(x[0],x[1]);
                }
            }

            // Finally, fill the lookup table:
            for (unsigned i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (pop.mcounts[i])
                        pop.mut_lookup.insert(pop.mutations[i].pos);
                }
        }

        /*!
          \brief Overload for metapopulation simulations
         */
        template <typename streamtype, typename sugarpop_t, typename reader_t,
                  typename diploid_reader_t = diploidIOplaceholder>
        inline typename std::
            enable_if<std::is_same<typename sugarpop_t::popmodel_t,
                                   sugar::METAPOP_TAG>::value,
                      result_type>::type
            operator()(sugarpop_t &pop, streamtype &buffer, const reader_t &rt,
                       const diploid_reader_t &dr = diploid_reader_t()) const
        {
            pop.clear();
            fwdpp_internal::scalar_reader reader;
            // Step 0: read N
            uint_t numNs;
            reader(buffer, &numNs);
            pop.Ns.resize(numNs);
            reader(buffer, &pop.Ns[0], numNs);
            // Step 1: write the mutations, diploids, gametes to the stream
            KTfwd::read_binary_metapop(pop.gametes, pop.mutations,
                                       pop.diploids, rt, buffer, dr);
            // update the mutation counts
            fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                            pop.mcounts);
            uint_t temp;
            reader(buffer, &temp);
            for (uint_t m = 0; m < temp; ++m)
                {
                    typename reader_t::result_type mm = rt(buffer);
                    pop.fixations.emplace_back(std::move(mm));
                }
            pop.fixation_times.resize(temp);
            if (temp)
                {
                    reader(buffer, &pop.fixation_times[0], temp);
                }

            // Finally, fill the lookup table:
            for (unsigned i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (pop.mcounts[i])
                        pop.mut_lookup.insert(pop.mutations[i].pos);
                }
        }
    };

    /*!
      Write a population to a gzFile in a binary format
    */
    struct gzserialize
    {
        using result_type = std::result_of<decltype (&gzwrite)(
            gzFile, void *, unsigned)>::type;
        /*
          \brief Call operator
          \note gzout must already be opened, and with a mode involving 'b'
         */
        template <typename sugarpop_t, typename writer_t,
                  typename diploid_writer_t = diploidIOplaceholder>
        inline result_type
        operator()(gzFile gzout, const sugarpop_t &pop, const writer_t &wt,
                   const diploid_writer_t &dw = diploid_writer_t()) const
        {
            std::ostringstream buffer;
            serialize()(buffer, pop, wt, dw);
            return gzwrite(gzout, buffer.str().c_str(),
                           unsigned(buffer.str().size()));
        }
    };

    /*!
      Read a population from a gzFile in binary format
      \deprecated
     */
    using gzdeserialize = deserialize;
}

#endif
