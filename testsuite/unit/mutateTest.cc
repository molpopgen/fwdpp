/*!
  \file mutateTest.cc
  \ingroup unit
  \brief Tests fwdpp::mutation, fwdpp::haploid_genome
*/

#include <config.h>
#include <fwdpp/forward_types.hpp>
#include <boost/test/unit_test.hpp>

// trivial ways to play with the fwdpp::mutation type
using mut = fwdpp::mutation;
using gtype = fwdpp::haploid_genome;

BOOST_AUTO_TEST_CASE(make_mutation_1)
{
    // Mutation at position 0.1, selection coefficient of 0
    mut m(0.1, 0.);

    BOOST_REQUIRE_EQUAL(m.pos, 0.1);
    BOOST_REQUIRE_EQUAL(m.neutral, true);
}

BOOST_AUTO_TEST_CASE(make_mutation_2)
{
    // Mutation at position 0.1, selection coefficient of 0,
    // dominance of 0.25
    mut m(0.1, 0, 0.25);

    BOOST_REQUIRE_EQUAL(m.pos, 0.1);
    BOOST_REQUIRE_EQUAL(m.neutral, true);
    BOOST_REQUIRE_EQUAL(m.h, 0.25);
}

BOOST_AUTO_TEST_CASE(copy_construct_1)
{
    mut m(0.1, 0., 1);
    mut m2(m);

    BOOST_REQUIRE(m == m2);
}

BOOST_AUTO_TEST_CASE(assign_1)
{
    mut m(0.1, 0., 1);
    mut m2 = m;

    BOOST_REQUIRE(m == m2);
}

/*
  This is a more realistic test:

  1. Add multiple mutations
  2. The positions of the mutations will not be nicely in order
*/
// BOOST_AUTO_TEST_CASE( add_N_mutations_1 )
// {
//   gtype g(1);

//   std::vector<double> next_mut_pos = { 2., 0.24, 0.25, 3, 1.999 };
//   decltype(next_mut_pos)::size_type i = 0;

//   std::vector<mut> mvector;
//   std::vector<fwdpp::uint_t> mcounts;
//   auto mut_recycling_bin = fwdpp::fwdpp_internal::make_mut_queue(mcounts);
//   //This is the mutation model.  Return a mutation at the next position
//   //Positions are deterministic for the sake of testing.
//   auto mmodel = [&next_mut_pos,&i]( decltype(mut_recycling_bin) &
//   rbin,std::vector<mut> & __mvector )
//     {
//       //mutations are all neutral
//       return
//       fwdpp::fwdpp_internal::recycle_mutation_helper(rbin,__mvector,next_mut_pos[i++],
//       0. );
//     };

//   fwdpp::fwdpp_internal::add_N_mutations_recycle(mut_recycling_bin,
// 						 mmodel,
// 						 next_mut_pos.size(),
// 						 mvector,
// 						 g);

//   BOOST_CHECK_EQUAL( mvector.size(), next_mut_pos.size() );
//   //neutral mutations should contain 5 things
//   BOOST_CHECK_EQUAL( g.mutations.size(), next_mut_pos.size() );
//   //selected mutations should be empty
//   BOOST_CHECK( g.smutations.empty() );
//   /*
//     This is the important test:

//     The library assumes that the iterators to mutations
//     stored by haploid_genomes are sorted w.r.to position.

//     The mutation functions take care of that.
//    */
//   BOOST_CHECK( std::is_sorted( g.mutations.cbegin(),
// 			       g.mutations.cend(),
// 			       [&mvector]( std::size_t i, std::size_t j ) {
// 				 return mvector[i].pos < mvector[j].pos;
// 			       } )
// 	       );

//   //We should also be able to find mutations by position
//   for( auto p : next_mut_pos )
//     {
//       BOOST_CHECK( std::find_if( g.mutations.begin(),
// 				 g.mutations.end(),
// 				 [&p,&mvector]( std::size_t i ) {
// 				   return mvector[i].pos == p;
// 				 } ) != g.mutations.end()
// 		   );
//     }
// }
