#ifndef __FWDPP_IO_HPP__
#define __FWDPP_IO_HPP__

#include <utility>
#include <zlib.h>

namespace KTfwd
{
  /*! \brief Write the population to a compact binary-format output file.
    Write the population to a compact binary-format output file.
    
    \param gametes  The vector of gametes
    \param mutations The list of mutations
    \param mw A function object taking a mutation and an ostreamtype as arguments. Must be prov
    \param buffer An ouptut stream into which the population is written.  This is the "return value" of the function.  The stream must support a write() function
    akin to those found in the std::ostream classes.
    
    \note If is often useful for buffer to be of type std::ostringstream to allow writing of the buffered data to C-style file handles/pointers, 
    in turn allowing file locking which speeds up performance on distributed file systems.
    
    \example diploid_binaryIO.cc
  */
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename mutation_type,
	    typename list_type_allocator,
	    template<typename,typename> class list_type,
	    typename mutation_writer,
	    typename ostreamtype>
  void write_binary_pop ( const vector_type< gamete_type, vector_type_allocator > * gametes,
			  list_type< mutation_type, list_type_allocator > * mutations,
			  const mutation_writer & mw,
			  ostreamtype & buffer);

  /*! \brief Write the metapopulation to a compact binary-format output file.
    Write the metapopulation to a compact binary-format output file.

    \param gametes  The vector of vector gametes for each deme in metapop
    \param mutations The list of mutations
    \param mw A function object taking a mutation and an ostreamtype as arguments. Must be prov
    \param buffer An ouptut stream into which the population is written.  This is the "return value" of the function.  The stream must support a write() function
    akin to those found in the std::ostream classes.

    \note If is often useful for buffer to be of type std::ostringstream to allow writing of the buffered data to C-style file handles/pointers, 
    in turn allowing file locking which speeds up performance on distributed file systems.
    \example diploid_binaryIO.cc
   */
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename vector_type_allocator2,
	    template<typename,typename> class vector_type2,
	    typename mutation_type,
	    typename list_type_allocator,
	    template<typename,typename> class list_type,
	    typename mutation_writer_type,
	    typename ostreamtype>
  void write_binary_metapop ( const vector_type2<vector_type< gamete_type, vector_type_allocator >, vector_type_allocator2 > * gametes,
			      list_type< mutation_type, list_type_allocator > * mutations,
			      const mutation_writer_type & mw,
			      ostreamtype & buffer);

  /*! \brief Read the population back from a binary-format file
    Read the population back from a binary-format file

    \param gametes Destination for the gametes
    \param mutations Destination for the mutations
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param in Input stream.  Must support .read() in a manner similar to std::istream types.

    \example diploid_binaryIO.cc
   */
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename mutation_type,
	    typename list_type_allocator,
	    template<typename,typename> class list_type,
	    typename mutation_reader,
	    typename istreamtype>
  void read_binary_pop ( vector_type< gamete_type, vector_type_allocator > * gametes,
			 list_type< mutation_type, list_type_allocator > * mutations,
			 const mutation_reader & mr,
			 istreamtype & in);

  /*! \brief Read the population back from a compressed binary-format file
    Read the population back from a compressed binary-format file

    \param gametes Destination for the gametes
    \param mutations Destination for the mutations
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param gzin An opend gzFile from zlib.

    \example diploid_gzbinaryIO.cc
   */
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_reader_type>
  void read_binary_pop (  gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			  mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
							  typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					       vector_type_allocator > * diploids,
			  const mutation_reader_type & mr,
			  gzFile gzin );

  /*! \brief Read the metapopulation back from a binary-format file
    Read the metapopulation back from a binary-format file

    \param gametes Destination for the gametes
    \param mutations Destination for the mutations
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param in Input stream.  Must support .read() in a manner similar to std::istream types.
   */
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename vector_type_allocator2,
	    template<typename,typename> class vector_type2,
	    typename mutation_type,
	    typename list_type_allocator,
	    template<typename,typename> class list_type,
	    typename mutation_reader,
	    typename istreamtype>
  void read_binary_metapop ( vector_type2< vector_type< gamete_type, vector_type_allocator >, vector_type_allocator2 > * gametes,
			     list_type< mutation_type, list_type_allocator > * mutations,
			     const mutation_reader & mr,
			     istreamtype & in);

  /*! \brief Read the metapopulation back from a compressed binary-format file
    Read the metapopulation back from a compressed binary-format file

    \param gametes Destination for the gametes
    \param mutations Destination for the mutations
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param gzin A gzFile (from zlib) opened for reading.
   */
  template< typename gamete_type,
	    typename vector_type_allocator,
	    template<typename,typename> class vector_type,
	    typename vector_type_allocator2,
	    template<typename,typename> class vector_type2,
	    typename mutation_type,
	    typename list_type_allocator,
	    template<typename,typename> class list_type,
	    typename mutation_reader>
  void read_binary_metapop ( vector_type2< vector_type< gamete_type, vector_type_allocator >, vector_type_allocator2 > * gametes,
			     list_type< mutation_type, list_type_allocator > * mutations,
			     const mutation_reader & mr,
			     gzFile gzin );

  /*! \brief Write population to binary-format file for individual-based simulations
    \param gametes  The vector of lists gametes for each deme in metapop
    \param mutations The list of mutations
    \param diploids The vector of diploids
    \param mw A function object taking a mutation and an ostreamtype as arguments. Must be prov
    \param buffer An ouptut stream into which the population is written.  This is the "return value" of the function.  The stream must support a write() function
    akin to those found in the std::ostream classes.

    \note If is often useful for buffer to be of type std::ostringstream to allow writing of the buffered data to C-style file handles/pointers, 
    in turn allowing file locking which speeds up performance on distributed file systems.
   */
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_writer_type,
	    typename ostreamtype>
  void write_binary_pop ( const gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			  mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  const diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
								typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
			  vector_type_allocator > * diploids,
			  const mutation_writer_type & mw,
			  ostreamtype & buffer);

  /*! \brief Read the population back from a binary-format file for individual-based simulations

    \param gametes Destination for the gametes
    \param mutations Destination for the mutations
    \param diploids Destination for the diploids
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param in Input stream.  Must support .read() in a manner similar to std::istream types.
   */
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_reader_type,
	    typename istreamtype>
  void read_binary_pop (  gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			  mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
							  typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					       vector_type_allocator > * diploids,
			  const mutation_reader_type & mr,
			  istreamtype & in);

  /*! \brief Read the population back from a compressed binary-format file for individual-based simulations
    Read the population back from a compressed binary-format file for individual-based simulations
    \param gametes Destination for the gametes
    \param mutations Destination for the mutations
    \param diploids Destination for the diploids
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param gzin A gzFile (from zlib) opened for reading.

    \example diploid_gzbinaryIO_ind.cc
   */
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_reader_type>
  void read_binary_pop (  gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			  mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
							  typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
					       vector_type_allocator > * diploids,
			  const mutation_reader_type & mr,
			  gzFile gzin );

  /*! \brief Write the metapopulation to a compact binary-format output file for individual-based simulations.
    Write the metapopulation to a compact binary-format output file.
    
    \param metapop  The vector of lists gametes for each deme in metapop
    \param mutations The list of mutations
    \param diploids The vector of vectors of diploids.
    \param mw A function object taking a mutation and an ostreamtype as arguments. Must be prov
    \param buffer An ouptut stream into which the population is written.  This is the "return value" of the function.  The stream must support a write() function
    akin to those found in the std::ostream classes.

    \note If is often useful for buffer to be of type std::ostringstream to allow writing of the buffered data to C-style file handles/pointers, 
    in turn allowing file locking which speeds up performance on distributed file systems.
    \example diploid_binaryIO_ind.cc
   */
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename metapop_vector_type_allocator,
	    template<typename,typename> class metapop_vector_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_writer_type,
	    typename ostreamtype>
  void write_binary_metapop ( const metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
			      metapop_vector_type_allocator> * metapop,
			      mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			      const diploid_vv_type < diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
										      typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
									   vector_type_allocator >,
			      diploid_vv_type_allocator > * diploids,
			      const mutation_writer_type & mw,
			      ostreamtype & buffer);

  /*! \brief Read the metapopulation back from a binary-format file for individual-based simulations
    Read the metapopulation back from a binary-format file

    \param metapop Destination for the gametes
    \param mutations Destination for the mutations
    \param diploids Destination for the diploids
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param in Input stream.  Must support .read() in a manner similar to std::istream types.
   */
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename metapop_vector_type_allocator,
	    template<typename,typename> class metapop_vector_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_reader_type,
	    typename istreamtype>
  void read_binary_metapop ( metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
			     metapop_vector_type_allocator> * metapop,
			     mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			     diploid_vv_type < diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
									       typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
								    vector_type_allocator >,
			     diploid_vv_type_allocator > * diploids,
			     const mutation_reader_type & mr,
			     istreamtype & in);

  /*! \brief Read the metapopulation back from a compressed binary-format file for individual-based simulations
    Read the metapopulation back from a compressed binary-format file

    \param metapop Destination for the gametes
    \param mutations Destination for the mutations
    \param diploids Destination for the diploids
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param gzin An opened gzFile from zlib.
   */
  template< typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename metapop_vector_type_allocator,
	    template<typename,typename> class metapop_vector_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_reader_type>
  void read_binary_metapop ( metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
			     metapop_vector_type_allocator> * metapop,
			     mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			     diploid_vv_type < diploid_vector_type< std::pair< typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator,
									       typename gamete_list_type< gamete_type, gamete_list_type_allocator >::iterator >,
								    vector_type_allocator >,
			     diploid_vv_type_allocator > * diploids,
			     const mutation_reader_type & mr,
			     gzFile gzin);

}
#endif
#include <fwdpp/IO.tcc>
