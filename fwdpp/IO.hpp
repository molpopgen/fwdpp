#ifndef __FWDPP_IO_HPP__
#define __FWDPP_IO_HPP__

#include <utility>
#include <fwdpp/io/scalar_serialization.hpp>

namespace fwdpp
{
    template <typename T> struct serialize_diploid
    {
        io::scalar_writer writer;
        serialize_diploid() : writer{} {}
        template <typename ostreamtype>
        inline void
        operator()(const T &dip, ostreamtype &buffer) const
        {
            writer(buffer, &dip.first);
            writer(buffer, &dip.second);
        }
    };

    template <typename T> struct deserialize_diploid
    {
        io::scalar_reader reader;
        deserialize_diploid() : reader{} {}
        template <typename istreamtype>
        inline void
        operator()(T &dip, istreamtype &buffer) const
        {
            typename T::first_type c;
            reader(buffer, &c);
            dip.first = c;
            reader(buffer, &c);
            dip.second = c;
        }
    };

    template <typename T> struct serialize_mutation
    {
        template <typename istreamtype>
        inline void
        operator()(const T &, istreamtype &) const
        {
            throw std::runtime_error(
                "serializtion not implemented for this mutation type");
        }
    };

    template <typename T> struct deserialize_mutation
    {
        template <typename ostreamtype>
        inline T
        operator()(const T &, ostreamtype &) const
        {
            throw std::runtime_error(
                "deserializtion not implemented for this mutation type");
        }
    };

    /*! \brief Write population to binary-format file for individual-based
      simulations
      \param gametes  The gamete container.
      \param mutations The container of mutations.
      \param diploids The vector of diploids.
      \param mw A function object taking a mutation and an ostreamtype as
      arguments. Must be provided by the library user.
      \param buffer An ouptut stream into which the population is written.
      This is the "return value" of the function.  The stream must either
      support a write() function
      akin to those found in the std::ostream classes or be a gzFile opened in
      'wb' mode.

      \note Using std::stringstream or std::ostringstream allows in-memory
      serialization for, e.g. MPI.
     */
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename ostreamtype>
    void write_binary_pop(const gcont_t &gametes, const mcont_t &mutations,
                          const dipvector_t &diploids, ostreamtype &buffer);

    /*! \brief Read the population back from a binary-format file for
      individual-based simulations

      \param gametes Destination for the gametes
      \param mutations Destination for the mutations
      \param diploids Destination for the diploids
      \param mr A function object to read in the mutation information. Takes an
      istreamtype as argument. Must be provided by library user.
      \param in Input stream. Must either support .read() in a manner similar
      to std::istream types or be a gzFile from zlib.
     */
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename istreamtype>
    void read_binary_pop(gcont_t &gametes, mcont_t &mutations,
                         dipvector_t &diploids, istreamtype &in);

    /*! \brief Write the population to a binary-format file for
      individual-based multilocus simulations.
      \param mlocus_gametes A container of gametes for a multilocus simulation
      \param mutations A container of mutation objects
      \param diploids A container of individuals in the simulation
      \param mw A function object taking a mutation and an ostreamtype as
      arguments. Must be provided by the library user.
      \param buffer An object whose public interface is compatible with
      std::ostream or is a gzFile
     */
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename ostreamtype>
    void write_binary_pop_mloc(const gcont_t &mlocus_gametes,
                               const mcont_t &mutations,
                               const dipvector_t &diploids,
                               ostreamtype &buffer);

    /*! \brief Read the population back from a binary-format file for
      individual-based multilocus simulations
      \param mlocus_gametes A container of gametes for a multilocus simulation
      \param mutations A mutation container
      \param diploids A container of individuals in the simulation
      \param mr A function object taking a input stream as argument, and reads
      a mutation object from the stream. Must be provided by the library user.
      \param in An object whose public interface is compatible with
      std::ostream or is a gzFile.
     */
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename istreamtype>
    void read_binary_pop_mloc(gcont_t &mlocus_gametes, mcont_t &mutations,
                              dipvector_t &diploids, istreamtype &in);

    /*! \brief Write the metapopulation to a compact binary-format output file
      for individual-based simulations.
      Write the metapopulation to a compact binary-format output file.

      \param gametes The gamete container.
      \param mutations THe mutation container.
      \param diploids The vector of vectors of diploids.
      \param mw A function object taking a mutation and an ostreamtype as
      arguments. Must be provided by the library user.
      \param buffer An ouptut stream into which the population is written.
      This is the "return value" of the function.  The stream must either
      support a write() function
      akin to those found in the std::ostream classes or be a gzFile opened in
      'wb' mode.

      \note Using std::stringstream or std::istringstream allows in-memory
      serialization for, e.g., MPI.
      \example diploid_binaryIO_ind.cc
     */
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename mutation_writer_type, typename ostreamtype>
    void write_binary_metapop(const gcont_t &gametes, const mcont_t &mutations,
                              const dipvector_t &diploids,
                              const mutation_writer_type &mw,
                              ostreamtype &buffer);

    /*! \brief Read the metapopulation back from a binary-format file for
      individual-based simulations
      Read the metapopulation back from a binary-format file
      \param metapop Destination for the gametes
      \param mutations Destination for the mutations
      \param diploids Destination for the diploids
      \param mr A function object to read in the mutation information. Takes an
      istreamtype as argument. Must be provided by library user.
      \param in Input stream.  Must either support .read() in a manner similar
      to std::istream types or be a gzFile from zlib.
     */
    template <typename gcont_t, typename mcont_t, typename dipvector_t,
              typename istreamtype>
    void read_binary_metapop(gcont_t &metapop, mcont_t &mutations,
                             dipvector_t &diploids, istreamtype &in);
}
#endif
#include <fwdpp/IO.tcc>
