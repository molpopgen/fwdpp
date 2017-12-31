#ifndef FWDPP_IO_SCALAR_SERIALIZATION_HPP__
#define FWDPP_IO_SCALAR_SERIALIZATION_HPP__

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace fwdpp
{
    namespace io
    {
        struct scalar_reader
        /// Serialize POD and contiguous arrays of POD
        {
            template <typename streamtype, typename T>
            inline void
            operator()(streamtype &i, T *__t, std::size_t n = 1) const
            {
                /*! \brief Read binary data
                 */
                i.read(reinterpret_cast<char *>(__t), n * sizeof(T));
            }
        };

        struct scalar_writer
        /// Deserialize POD and contiguous arrays of POD
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
        };
    }
}

#endif
