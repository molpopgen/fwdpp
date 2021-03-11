/*! \file GSLrng_t.hpp
  \brief Wrapper for gsl_rng *
 */
#ifndef __FWDPP_SUGAR_GSLRNG_T_HPP__
#define __FWDPP_SUGAR_GSLRNG_T_HPP__

#include <stdexcept>
//#include <cstdio>// FOR SOME POSIX ONLY FUNCTIONS!
#include <fwdpp/gsl/tags.hpp>
#include <fwdpp/gsl/deleter.hpp>

namespace fwdpp
{

    //! Distpatch tag to signal GSLrng_t to instantiate in terms of
    //! gsl_rng_mt19937
    using GSL_RNG_MT19937
        = gsl::GSL_RNG_TYPE_TAG<gsl::GSL_RNG_TYPE::MT19937>;
    //! Distpatch tag to signal GSLrng_t to instantiate in terms of
    //! gsl_rng_taus2
    using GSL_RNG_TAUS2 = gsl::GSL_RNG_TYPE_TAG<gsl::GSL_RNG_TYPE::TAUS2>;

    /*!
      \brief A wrapper around gsl_rng * objects.

      The template instantiation type must be a model of
      fwdpp::gsl::GSL_RNG_TYPE_TAG, which specifies the
      gsl_rng type.
      This type holds an object of type fwdpp::gsl::gsl_rng_ptr_t,
      which is a smart pointer that manages freeing the gsl_rng * upon
      destruction.
     */
    template <typename T> class GSLrng_t
    {
      private:
        gsl::gsl_rng_ptr_t setup(GSL_RNG_MT19937)
        {
            return gsl::gsl_rng_ptr_t(gsl_rng_alloc(gsl_rng_mt19937),
                                        [](gsl_rng *r) { gsl_rng_free(r); });
        }

        gsl::gsl_rng_ptr_t setup(GSL_RNG_TAUS2)
        {
            return gsl::gsl_rng_ptr_t(gsl_rng_alloc(gsl_rng_taus2),
                                        [](gsl_rng *r) { gsl_rng_free(r); });
        }

        gsl::gsl_rng_ptr_t
        setup(const gsl_rng *r)
        {
            return gsl::gsl_rng_ptr_t(gsl_rng_clone(r),
                                        [](gsl_rng *r) { gsl_rng_free(r); });
        }

        //! Smart pointer wrapping the gsl_rng *
        gsl::gsl_rng_ptr_t r;

      public:
        //! Typedef for RNG type, if needed
        using rngtype = T;

        //! Construct with a seed
        GSLrng_t(const unsigned long seed) noexcept : r(setup(T()))
        {
            gsl_rng_set(r.get(), seed);
        }

        ~GSLrng_t() = default;
        GSLrng_t(const GSLrng_t &rng) = delete;
		/// Move constructor
        GSLrng_t(GSLrng_t &&) noexcept = default;
        GSLrng_t &operator=(GSLrng_t &) = delete;
		/// Move assignment
        GSLrng_t &operator=(GSLrng_t &&) noexcept = default;

        //! Return underlying pointer
        const gsl_rng *
        get() const
        {
            return r.get();
        }

        // /*!
        //   Write gsl_rng state to memory buffer.

        //   Throws std::runtime_error if there are problems
        // */
        // template<typename ostream_t>
        // void serialize(ostream_t & o) const
        // {
        //   /*
        // 	We use C-style NULL in lieu of nullptr
        // 	b/c we cannot be sure what the POSIX functions
        // 	are doing...
        //   */
        //   FILE * stream = NULL;
        //   char * buf = NULL;
        //   std::size_t len;
        //   stream = open_memstream(&buf,&len);
        //   auto rv = gsl_rng_fwrite(stream,this->get());
        //   if( rv ) //gsl_rng_fwrite returns 0 upon success
        // 	{
        // 	  fclose(stream);
        // 	  if(buf != NULL) free(buf);
        // 	  throw std::runtime_error("error calling gsl_rng_fwrite");
        // 	}
        //   fclose(stream);
        //   o.write( reinterpret_cast<char*>(&len),sizeof(std::size_t) );
        //   o.write( buf, len*sizeof(char) );
        //   free(buf);
        // }

        // //! Read RNG state from a memory buffer
        // template<typename istream_t>
        // void deserialize(istream_t & i)
        // {
        //   std::size_t len;
        //   i.read(reinterpret_cast<char*>(&len),sizeof(std::size_t));
        //   std::unique_ptr<char> buffer(new char[len]);
        //   i.read(buffer.get(),len*sizeof(char));
        //   FILE * stream = fmemopen(buffer.get(),len*sizeof(char),"r");
        //   auto rv = gsl_rng_fread(stream,r.get());
        //   fclose(stream);
        // }
    };

	/// \typedef GSLrng_mt
	/// Typedef for mersenne twister
	/// \version 0.7.4 Added to fwdpp
	using GSLrng_mt = GSLrng_t<GSL_RNG_MT19937>;
}

#endif
