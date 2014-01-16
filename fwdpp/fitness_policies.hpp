#ifndef _KTFWD_FITNESS_POLICIES_HPP_
#define _KTFWD_FITNESS_POLICIES_HPP_

namespace KTfwd
{
  /*! \brief What to do with a diploid that is homozygous for a specific mutation under multiplicative fitness
    \param fitness  The diploid's current fitness
    \param g1 Pointer to one parent's mutation for which the diploid is homozygous
    \param scaling  See documentation of return value.  Allows you to simulate homozygote fitnesses of 1+s or 1+2s, etc.
    \return fitness *= (1 + scaling*(g1->s)); 
   */
  struct multiplicative_fitness_updater_hom
  {
    typedef void result_type;
    template<typename iterator_type>
    inline void operator()(double & fitness, const iterator_type & m1,const double & scaling = 2.) const
    {
      fitness *= ( 1. + scaling*m1->s );
    }
  };
  
  /*! \brief What to do with a diploid that is heterozygous for a specific mutation under multiplicative fitness
    \param fitness  The diploid's current fitness
    \param g1 Pointer to one parent's mutation for which the diploid is homozygous
    \param scaling  See documentation of return value
    \return fitness *= (1 + scaling*(g1->s)*(g1->h)); 
    \note The value_type of m1 must have data types s and h!!
   */
  struct multiplicative_fitness_updater_het
  {
    typedef void result_type;
    template<typename iterator_type>
    inline void operator()(double & fitness, const iterator_type & m1) const
    {
      fitness *= ( 1. + m1->s*m1->h );
    }
  };

  /*! \brief What to do with a diploid that is homozygous for a specific mutation under additive fitness
    \param fitness  The diploid's current fitness
    \param g1 Pointer to one parent's mutation for which the diploid is homozygous
    \param scaling  See documentation of return value.  Allows you to simulate homozygote fitnesses of 1+s or 1+2s, etc.
    \return fitness += (1 + scaling*(g1->s)); 
   */
  struct additive_fitness_updater_hom
  {
    typedef void result_type;
    template<typename iterator_type>
    inline void operator()(double & fitness, const iterator_type & m1,const double & scaling = 2.) const
    {
      fitness += ( scaling*m1->s );
    }
  };
  
  /*! \brief What to do with a diploid that is heterozygous for a specific mutation under additive fitness
    \param fitness  The diploid's current fitness
    \param g1 Pointer to one parent's mutation for which the diploid is homozygous
    \param scaling  See documentation of return value
    \return fitness += (1 + scaling*(g1->s)*(g1->h)); 
    \note The value_type of m1 must have data types s and h!!
   */
  struct additive_fitness_updater_het
  {
    typedef void result_type;
    template<typename iterator_type>
    inline void operator()(double & fitness, const iterator_type & m1) const
    {
      fitness += ( m1->s*m1->h );
    }
  };
}
#endif
