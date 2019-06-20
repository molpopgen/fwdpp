#ifndef __FWDPP_INTERNAL_GSL_DISCRETE_HPP__
#define __FWDPP_INTERNAL_GSL_DISCRETE_HPP__

#include <gsl/gsl_randist.h>
#include <memory>

namespace fwdpp
{
    namespace detail
    {
        /*!
          Custom deleter for std::uniq_ptr
		*/
        struct gsl_ran_discrete_t_deleter
        {
            void
            operator()(gsl_ran_discrete_t *l) noexcept
            {
                gsl_ran_discrete_free(l);
            }
        };
    } // namespace detail

    /*!
	  \warning Can only be put into a vector by push_back/emplace back
	  because of constraints on std::unique_ptr assignment
	 */
    using gsl_ran_discrete_t_ptr
        = std::unique_ptr<gsl_ran_discrete_t, detail::gsl_ran_discrete_t_deleter>;
} // namespace fwdpp

#endif
