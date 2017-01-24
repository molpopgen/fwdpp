#ifndef __FWDPP_INTERNAL_IOHELP_HPP__
#define __FWDPP_INTERNAL_IOHELP_HPP__

/*
  Mechanics of data serialization
  The various write_binary_pop and read_binary pop
  functions rely on these implementations
*/
#include <zlib.h>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <fwdpp/forward_types.hpp>

namespace KTfwd
{
    namespace fwdpp_internal
    {

        struct scalar_reader
        {
            template <typename streamtype, typename T>
            inline void
            operator()(streamtype &i, T *__t, std::size_t n = 1) const
            {
                /*! \brief Read binary data
                 */
                i.read(reinterpret_cast<char *>(__t), n * sizeof(T));
            }
            template <typename T>
            inline void
            operator()(gzFile &gzin, T *__t, std::size_t n = 1) const
            {
                /*! \brief Read binary data
                 */
                gzread(gzin, __t, n * sizeof(T));
            }
        };

        struct scalar_writer
        {
            using result_type = std::uint64_t;
            template <typename streamtype, typename T>
            inline result_type
            operator()(streamtype &i, T *__t, std::size_t n = 1) const
            {
                /*! \brief Write binary data
                 * \throw std::runtime_error
                 */
                if (!i)
                    {
                        throw std::runtime_error("serialization error on line "
                                                 + std::to_string(__LINE__)
                                                 + " of "
                                                 + std::string(__FILE__));
                    }
                i.write(reinterpret_cast<const char *>(__t), n * sizeof(T));
                if (!i)
                    {
                        throw std::runtime_error("serialization error on line "
                                                 + std::to_string(__LINE__)
                                                 + " of "
                                                 + std::string(__FILE__));
                    }
                return result_type(n * sizeof(T));
            }
            template <typename T>
            inline result_type
            operator()(gzFile &gzout, T *__t, std::size_t n = 1) const
            {
                /*! \brief Write binary data
                 * \throw std::runtime_error
                 */
                auto rv = gzwrite(gzout, __t, n * sizeof(T));
                if (!rv)
                    {
                        throw std::runtime_error("serialization error on line "
                                                 + std::to_string(__LINE__)
                                                 + " of "
                                                 + std::string(__FILE__));
                    }
                return result_type(rv);
            }
        };

        struct write_mutations
        {
            template <typename mutation_type,
                      typename container_type_allocator,
                      template <typename, typename> class container_type,
                      typename mutation_writer_type, typename ostreamtype>
            void
            operator()(
                const container_type<mutation_type, container_type_allocator>
                    &mutations,
                const mutation_writer_type &mw, ostreamtype &buffer) const
            {
                std::size_t MUTNO = mutations.size();
                scalar_writer()(buffer, &MUTNO);
                // write the mutation data to the buffer
                for (const auto &m : mutations)
                    mw(m, buffer);
            }
        };

        struct write_haplotypes
        {
            template <typename gamete_type, typename... gamete_cont_t_details,
                      template <typename, typename...> class gamete_cont_t,
                      typename ostreamtype>
            void
            operator()(const gamete_cont_t<gamete_type,
                                           gamete_cont_t_details...> &gametes,
                       ostreamtype &buffer) const
            {
                std::size_t N = gametes.size();
                scalar_writer writer;
                writer(buffer, &N);
                for (const auto &g : gametes)
                    {
                        writer(buffer, &g.n);
                        std::size_t nm = g.mutations.size();
                        writer(buffer, &nm);
                        if (nm)
                            {
                                writer(buffer, g.mutations.data(), nm);
                            }
                        nm = g.smutations.size();
                        writer(buffer, &nm);
                        if (nm)
                            {
                                writer(buffer, g.smutations.data(), nm);
                            }
                    }
            }
        };

        struct read_mutations
        {
            template <typename mutation_type,
                      typename container_type_allocator,
                      template <typename, typename> class container_type,
                      typename mutation_reader, typename istreamtype>
            void
            operator()(container_type<mutation_type, container_type_allocator>
                           &mutations,
                       const mutation_reader &mr, istreamtype &in) const
            {
                std::size_t NMUTS;
                scalar_reader()(in, &NMUTS);
                for (uint_t i = 0; i < NMUTS; ++i)
                    {
                        mutations.emplace_back(mr(in));
                    }
            }
        };

        struct read_haplotypes
        {
            template <typename gamete_type, typename container_type_allocator,
                      template <typename, typename> class container_type,
                      typename istreamtype>
            void
            operator()(
                container_type<gamete_type, container_type_allocator> &gametes,
                istreamtype &in) const
            {
                scalar_reader reader;
                std::size_t NHAPS;
                reader(in, &NHAPS);
                uint_t N;
                std::size_t nm;
                for (uint_t i = 0; i < NHAPS; ++i)
                    {
                        reader(in, &N);
                        gamete_type g(N);
                        reader(in, &nm);
                        if (nm)
                            {
                                g.mutations.resize(nm);
                                reader(in, g.mutations.data(), nm);
                            }
                        reader(in, &nm);
                        if (nm)
                            {
                                g.smutations.resize(nm);
                                reader(in, g.smutations.data(), nm);
                            }
                        gametes.emplace_back(std::move(g));
                    }
            }
        };
    }
}

#endif
