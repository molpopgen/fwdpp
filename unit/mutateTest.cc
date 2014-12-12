#define BOOST_TEST_MODULE mutateTest
#define BOOST_TEST_DYN_LINK 

#include <config.h>
#include <fwdpp/diploid.hh>
#include <boost/test/unit_test.hpp>
#include <unistd.h>
#include <iterator>
#include <functional>
#include <list>
//trivial ways to play with the KTfwd::mutation type
using mut = KTfwd::mutation;
using gtype = KTfwd::gamete;

BOOST_AUTO_TEST_CASE( make_mutation_1 )
{
  //Mutation at position 0.1, selection coefficient of 0, count in population of 1
  mut m(0.1,0.,1);

  BOOST_REQUIRE_EQUAL(m.pos,0.1);
  BOOST_REQUIRE_EQUAL(m.n,1);
  BOOST_REQUIRE_EQUAL(m.neutral,true);
}

BOOST_AUTO_TEST_CASE( make_mutation_2 )
{
  //Mutation at position 0.1, selection coefficient of 0, count in population of 1
  //dominance of 0.25
  mut m(0.1,0.,1,0.25);

  BOOST_REQUIRE_EQUAL(m.pos,0.1);
  BOOST_REQUIRE_EQUAL(m.n,1);
  BOOST_REQUIRE_EQUAL(m.neutral,true);
  BOOST_REQUIRE_EQUAL(m.h,0.25);
}

BOOST_AUTO_TEST_CASE( copy_construct_1 )
{
  mut m(0.1,0.,1);
  mut m2(m);

  BOOST_REQUIRE(m == m2);
}

BOOST_AUTO_TEST_CASE( assign_1 )
{
  mut m(0.1,0.,1);
  mut m2 = m;

  BOOST_REQUIRE(m == m2);
}

//Tests the internal machinery for adding a mutation to a gamete
BOOST_AUTO_TEST_CASE( add_mutation_1 )
{
  //create a gamete at frequency 1
  gtype g(1);

  //create a neutral mutation and add it to the "mutation pool"
  std::list<mut> mlist;
  auto mitr = mlist.insert(mlist.end(),mut(0.1,0.,1));

  KTfwd::fwdpp_internal::add_new_mutation(mitr,g);

  //So now, there is 1 neutral mutation and no selected mutations
  BOOST_CHECK_EQUAL( g.mutations.size(), 1 );
  BOOST_CHECK( g.smutations.empty() );

  //let's put in a selected mutation
  mitr = mlist.insert(mlist.end(),mut(0.1,-2.,1));
  KTfwd::fwdpp_internal::add_new_mutation(mitr,g);

  BOOST_CHECK_EQUAL( g.mutations.size(), 1 );
  BOOST_CHECK_EQUAL( g.smutations.size(), 1 );
}

/*
  This is a more realistic test:

  1. Add multiple mutations
  2. The positions of the mutations will not be nicely in order
*/
BOOST_AUTO_TEST_CASE( add_N_mutations_1 )
{
  gtype g(1);

  std::vector<double> next_mut_pos = { 2., 0.24, 0.25, 3, 1.999 };
  decltype(next_mut_pos)::size_type i = 0;

  std::list<mut> mlist;

  //This is the mutation model.  Return a mutation at the next position
  //Positions are deterministic for the sake of testing.
  auto mmodel = [&next_mut_pos,&i]( std::list<mut> * __mlist )
    {
      //mutations are all neutral
      return mut( next_mut_pos[i++], 0., 1 );
    };

  KTfwd::fwdpp_internal::add_N_mutations(mmodel,
					 /*
					   The policy below is equivalent to:
					   [](const mut & __mut, std::list<mut> * __mlist){ return __mlist->insert(__mlist->end(),__mut); }
					  */
					 std::bind(KTfwd::insert_at_end<mut,std::list<mut> >,std::placeholders::_1,std::placeholders::_2),
					 next_mut_pos.size(),
					 &mlist,
					 g);
  BOOST_CHECK_EQUAL( mlist.size(), next_mut_pos.size() );
  //neutral mutations should contain 5 things
  BOOST_CHECK_EQUAL( g.mutations.size(), next_mut_pos.size() );
  //selected mutations should be empty
  BOOST_CHECK( g.smutations.empty() );
  /*
    This is the important test:
    
    The library assumes that the iterators to mutations
    stored by gametes are sorted w.r.to position.

    The mutation functions take care of that.
   */
  BOOST_CHECK( std::is_sorted( g.mutations.cbegin(),
			       g.mutations.cend(),
			       []( std::list<mut>::iterator i,
				   std::list<mut>::iterator j ) {
				 return i->pos < j->pos;
			       } )
	       );

  //We should also be able to find mutations by position
  for( auto p : next_mut_pos )
    {
      BOOST_CHECK( std::find_if( g.mutations.begin(),
				 g.mutations.end(),
				 [&p]( std::list<mut>::iterator i ) {
				   return i->pos == p;
				 } ) != g.mutations.end()
		   );
    }
}

/*
  Same as previous test, but all mutations at the same position.
  This is allowed,  but that doesn't meant it is a good idea 
  unless you really know what you're doing
*/
BOOST_AUTO_TEST_CASE( add_N_mutations_2 )
{
  gtype g(1);

  std::vector<double> next_mut_pos = { 1.,1.,1.,1. };
  decltype(next_mut_pos)::size_type i = 0;

  std::list<mut> mlist;

  //This is the mutation model.  Return a mutation at the next position
  //Positions are deterministic for the sake of testing.
  auto mmodel = [&next_mut_pos,&i]( std::list<mut> * __mlist )
    {
      //mutations are all neutral
      return mut( next_mut_pos[i++], 0., 1 );
    };

  KTfwd::fwdpp_internal::add_N_mutations(mmodel,
					 /*
					   The policy below is equivalent to:
					   [](const mut & __mut, std::list<mut> * __mlist){ return __mlist->insert(__mlist->end(),__mut); }
					  */
					 std::bind(KTfwd::insert_at_end<mut,std::list<mut> >,std::placeholders::_1,std::placeholders::_2),
					 next_mut_pos.size(),
					 &mlist,
					 g);
  BOOST_CHECK_EQUAL( mlist.size(), next_mut_pos.size() );
  //neutral mutations should contain 5 things
  BOOST_CHECK_EQUAL( g.mutations.size(), next_mut_pos.size() );
  //selected mutations should be empty
  BOOST_CHECK( g.smutations.empty() );
  /*
    This is the important test:
    
    The library assumes that the iterators to mutations
    stored by gametes are sorted w.r.to position.

    The mutation functions take care of that.
   */
  BOOST_CHECK( std::is_sorted( g.mutations.cbegin(),
			       g.mutations.cend(),
			       []( std::list<mut>::iterator i,
				   std::list<mut>::iterator j ) {
				 return i->pos < j->pos;
			       } )
	       );
}
