#ifndef _INSERTION_POLICIES_HPP_
#define _INSERTION_POLICIES_HPP_

#include <algorithm>
#include <fwdpp/fwd_functional.hpp>
namespace KTfwd
{
    /*! \brief     An insertion policy
      Uses object forwarding to add an object of
      type T into container of type cT (e.g., container<T>)
     */
    struct emplace_back
    {
        /// \return the index where the emplaced object is located
        template <typename T, typename cT>
        inline std::size_t
        operator()(T &&t, cT &ct) const
        {
            ct.emplace_back(std::forward<T>(t));
            return ct.size() - 1;
        }
    };
}
#endif /* _INSERTION_POLICIES_HPP_ */
