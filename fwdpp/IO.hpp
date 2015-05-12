#ifndef __FWDPP_IO_HPP__
#define __FWDPP_IO_HPP__

#include <utility>

namespace KTfwd
{
  /*! \brief Write the population to a compact binary-format output file.
    Write the population to a compact binary-format output file.
    
    \param gametes  The vector of gametes
    \param mutations The list of mutations
    \param mw A function object taking a mutation and an ostreamtype as arguments. Must be provided by the library user.
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
			  const list_type< mutation_type, list_type_allocator > * mutations,
			  const mutation_writer & mw,
			  ostreamtype & buffer);

  /*! \brief Write the metapopulation to a compact binary-format output file.
    Write the metapopulation to a compact binary-format output file.

    \param gametes  The vector of vector gametes for each deme in metapop
    \param mutations The list of mutations
    \param mw A function object taking a mutation and an ostreamtype as arguments. Must be provided by the library user.
    \param buffer An ouptut stream into which the population is written.  This is the "return value" of the function.  The stream must support a write() function
    akin to those found in the std::ostream classes.

    \note If is often useful for buffer to be of type std::ostringstream to allow writing of the buffered data to C-style file handles/pointers, 
    in turn allowing file locking which speeds up performance on distributed file systems.
   */
  template<typename mutation_type,
	   typename list_type_allocator,
	   template<typename,typename> class list_type,
	   template<typename,typename> class gamete_type,
	   typename vector_type_allocator,
	   template<typename,typename> class vector_type,
	   typename vector_type_allocator2,
	   template<typename,typename> class vector_type2,
	   typename mutation_writer_type,
	   typename ostreamtype>
  void write_binary_metapop ( const vector_type2<vector_type< gamete_type<mutation_type,list_type<mutation_type,list_type_allocator> >, vector_type_allocator >, vector_type_allocator2 > * gametes,
			      const list_type< mutation_type, list_type_allocator > * mutations,
			      const mutation_writer_type & mw,
			      ostreamtype & buffer);
  
  /*! \brief Read the population back from a binary-format file
    Read the population back from a binary-format file

    \param gametes Destination for the gametes
    \param mutations Destination for the mutations
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param in Input stream.  Must either support .read() in a manner similar to std::istream types or be a gzFile from zlib.
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


  /*! \brief Read the metapopulation back from a binary-format file
    Read the metapopulation back from a binary-format file

    \param gametes Destination for the gametes
    \param mutations Destination for the mutations
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param in Input stream.  Must either support .read() in a manner similar to std::istream types or be a gzFile from zlib.
   */
   template<typename mutation_type,
	   typename list_type_allocator,
	   template<typename,typename> class list_type,
	   template<typename,typename> class gamete_type,
	   typename vector_type_allocator,
	   template<typename,typename> class vector_type,
	   typename vector_type_allocator2,
	   template<typename,typename> class vector_type2,
	   typename mutation_reader,
	   typename istreamtype>
  void read_binary_metapop ( vector_type2< vector_type< gamete_type<mutation_type,list_type<mutation_type,list_type_allocator> >, vector_type_allocator >, vector_type_allocator2 > * gametes,
			     list_type< mutation_type, list_type_allocator > * mutations,
			     const mutation_reader & mr,
			     istreamtype & in);

  /*!
    \brief placeholder IO policy for standard diploid types
  */
  struct diploidIOplaceholder
    {
      using result_type = void;
      /*!
	Does nothing
      */
      template< typename iterator, typename streamtype >
      inline result_type operator()( iterator , streamtype & ) const
      {
	//Does nothing!
      }
    };
  
  /*! \brief Write population to binary-format file for individual-based simulations
    \param gametes  The vector of lists gametes for each deme in metapop
    \param mutations The list of mutations
    \param diploids The vector of diploids
    \param mw A function object taking a mutation and an ostreamtype as arguments. Must be provided by the library user.
    \param buffer An ouptut stream into which the population is written.  This is the "return value" of the function.  The stream must support a write() function
    akin to those found in the std::ostream classes.

    \note If is often useful for buffer to be of type std::ostringstream to allow writing of the buffered data to C-style file handles/pointers, 
    in turn allowing file locking which speeds up performance on distributed file systems.
   */
  template< typename diploid_geno_t,
	    typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t = diploidIOplaceholder>
  void write_binary_pop ( const gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			  const mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  const diploid_vector_type< diploid_geno_t,vector_type_allocator > * diploids,
			  const mutation_writer_type & mw,
			  ostreamtype & buffer,
			  const diploid_writer_t & dw = diploid_writer_t());



  /*! \brief Read the population back from a binary-format file for individual-based simulations

    \param gametes Destination for the gametes
    \param mutations Destination for the mutations
    \param diploids Destination for the diploids
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param in Input stream. Must either support .read() in a manner similar to std::istream types or be a gzFile from zlib.
   */
  template< typename diploid_geno_t,
	    typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename mutation_reader_type,
	    typename istreamtype,
	    typename diploid_reader_t = diploidIOplaceholder>
  void read_binary_pop ( gamete_list_type< gamete_type, gamete_list_type_allocator > * gametes,
			 mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			 diploid_vector_type< diploid_geno_t,vector_type_allocator > * diploids,
			 const mutation_reader_type & mr,
			 istreamtype & in,
			 const diploid_reader_t & dr = diploid_reader_t());

  /*! \brief Write the population to a binary-format file for individual-based multilocus simulations.
    \param mlocus_gametes A container of gametes for a multilocus simulation
    \param mutations A linked list of mutation objects
    \param diploids A container of individuals in the simulation
    \param mw A function object taking a mutation and an ostreamtype as arguments. Must be provided by the library user.
    \param buffer An object whose public interface is compatible with std::ostream or is a gzFile
   */
  template< typename diploid_geno_t,
	    typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mlocus_vector_type_allocator,
	    template<typename,typename> class mlocus_vector_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_writer_type,
	    typename ostreamtype,
	    typename diploid_writer_t = diploidIOplaceholder>
  void write_binary_pop ( const mlocus_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >, mlocus_vector_type_allocator> * mlocus_gametes,
			  const mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			  const diploid_vv_type < diploid_vector_type< diploid_geno_t ,vector_type_allocator >,diploid_vv_type_allocator > * diploids,
			  const mutation_writer_type & mw,
			  ostreamtype & buffer,
			  const diploid_writer_t & dw = diploid_writer_t());

  /*! \brief Read the population back from a binary-format file for individual-based multilocus simulations
    \param mlocus_gametes A container of gametes for a multilocus simulation
    \param mutations A linked list of mutation objects
    \param diploids A container of individuals in the simulation
    \param mr A function object taking a input stream as argument, and reads a mutation object from the stream. Must be provided by the library user.
    \param in An object whose public interface is compatible with std::ostream or is a gzFile.
   */
  template< typename diploid_geno_t,
	    typename gamete_type,
	    typename gamete_list_type_allocator,
	    template<typename,typename> class gamete_list_type,
	    typename mlocus_vector_type_allocator,
	    template<typename,typename> class mlocus_vector_type,
	    typename mutation_type,
	    typename mutation_list_type_allocator,
	    template<typename,typename> class mutation_list_type,
	    typename vector_type_allocator,
	    template<typename,typename> class diploid_vector_type,
	    typename diploid_vv_type_allocator,
	    template<typename,typename> class diploid_vv_type,
	    typename mutation_reader_type,
	    typename istreamtype,
	    typename diploid_reader_t = diploidIOplaceholder>
  void read_binary_pop ( mlocus_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >, mlocus_vector_type_allocator> * mlocus_gametes,
			 mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			 diploid_vv_type < diploid_vector_type< diploid_geno_t , vector_type_allocator >, diploid_vv_type_allocator > * diploids,
			 const mutation_reader_type & mr,
			 istreamtype & in,
			 const diploid_reader_t & dr = diploid_reader_t());

  /*! \brief Write the metapopulation to a compact binary-format output file for individual-based simulations.
    Write the metapopulation to a compact binary-format output file.
    
    \param metapop  The vector of lists gametes for each deme in metapop
    \param mutations The list of mutations
    \param diploids The vector of vectors of diploids.
    \param mw A function object taking a mutation and an ostreamtype as arguments. Must be provided by the library user.
    \param buffer An ouptut stream into which the population is written.  This is the "return value" of the function.  The stream must support a write() function
    akin to those found in the std::ostream classes.

    \note If is often useful for buffer to be of type std::ostringstream to allow writing of the buffered data to C-style file handles/pointers, 
    in turn allowing file locking which speeds up performance on distributed file systems.
    \example diploid_binaryIO_ind.cc
   */
  template< typename diploid_geno_t,
	    typename gamete_type,
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
	    typename ostreamtype,
	    typename diploid_writer_t = diploidIOplaceholder>
  void write_binary_metapop ( const metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
			      metapop_vector_type_allocator> * metapop,
			      const mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			      const diploid_vv_type < diploid_vector_type< diploid_geno_t , vector_type_allocator >, diploid_vv_type_allocator > * diploids,
			      const mutation_writer_type & mw,
			      ostreamtype & buffer,
			      const diploid_writer_t & dw = diploid_writer_t());

  /*! \brief Read the metapopulation back from a binary-format file for individual-based simulations
    Read the metapopulation back from a binary-format file
    \param metapop Destination for the gametes
    \param mutations Destination for the mutations
    \param diploids Destination for the diploids
    \param mr A function object to read in the mutation information. Takes an istreamtype as argument. Must be provided by library user.
    \param in Input stream.  Must either support .read() in a manner similar to std::istream types or be a gzFile from zlib.
   */
  template< typename diploid_geno_t,
	    typename gamete_type,
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
  	    typename istreamtype,
	    typename diploid_reader_t = diploidIOplaceholder>
  void read_binary_metapop ( metapop_vector_type< gamete_list_type< gamete_type, gamete_list_type_allocator >,
			     metapop_vector_type_allocator> * metapop,
			     mutation_list_type< mutation_type, mutation_list_type_allocator > * mutations,
			     diploid_vv_type < diploid_vector_type<diploid_geno_t ,vector_type_allocator >, diploid_vv_type_allocator > * diploids,
			     const mutation_reader_type & mr,
			     istreamtype & in,
			     const diploid_reader_t & dr = diploid_reader_t());

}
#endif
#include <fwdpp/IO.tcc>
