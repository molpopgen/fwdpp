#ifndef FWDPP_INTERNAL_VOID_T_HPP__
#define FWDPP_INTERNAL_VOID_T_HPP__

namespace fwdpp
{
    namespace traits
    {
        namespace internal
        {
            // Based on
            // http://stackoverflow.com/questions/11813940/possible-to-use-type-traits-sfinae-to-find-if-a-class-defines-a-member-typec
            template <typename...> struct void_t
            {
                typedef void type;
            };
		}
	}
}
#endif
