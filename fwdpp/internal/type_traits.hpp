#ifndef FWDPP_INTERNAL_TYPE_TRAITS_HPP
#define FWDPP_INTERNAL_TYPE_TRAITS_HPP

namespace KTfwd
{
    namespace traits
    {
        namespace internal
        {
            // Based on
            // http://stackoverflow.com/questions/11813940/possible-to-use-type-traits-sfinae-to-find-if-a-class-defines-a-member-typec

            template <class T> struct void_t
            {
                typedef void type;
            };

/*
   This macro is used to generate
   types able to determine if a class type
   has a particular member type name
*/
#define HAS_TYPE(NAME)                                                        \
    template <typename, typename = void> struct has_##NAME : std::false_type  \
    {                                                                         \
    };                                                                        \
    template <typename T>                                                     \
    struct has_##NAME<T, typename void_t<typename T::NAME>::type>             \
        : std::true_type                                                      \
    {                                                                         \
    };

            HAS_TYPE(gamete_tag);
            HAS_TYPE(mutation_type);
            HAS_TYPE(first_type);
            HAS_TYPE(second_type);
        }
    }
}

#endif
