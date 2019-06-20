#include <iostream>
#include <config.h>
#include <cmath>
#include <queue>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <fwdpp/ts/table_simplifier.hpp>
#include <fwdpp/ts/generate_offspring.hpp>
#include <fwdpp/ts/get_parent_ids.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <boost/test/unit_test.hpp>


/* Let's describe the fixture setup, expected outputs, etc.:
 *
 * The crossover functions in the fixture will assign the following numbers
 * of xover events w/in and b/w loci:
 *
 * Locus w/in b/w*
 * 0     0    N/A
 * 1     1    0
 * 2     0    1
 * 3     1    2
 *
 * *The entry for b/w loci refers to the number of xover events
 * b/w locus i-1 and i.
 *
 * An xover w/in locus i occurs at position i + 0.5
 *
 * Let's draw the 4 loci in two parent chromosomes:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++++++ ++++++++ ++++++++   
 * -------- -------- -------- --------   
 *
 * The way multi-locus recombination works is that we apply
 * each recombination event from left to right.  Internally, this
 * is not literally what happens, but it is our mental image of 
 * what we're trying to accomplish
 *
 * We will consider the case of NO INITIAL SWAPPING of either
 * parent gamete, and "build" the first offspring gamete.
 *
 * There are no x-overs w/in locus 1, so after processing it,
 * we should have:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++++++ ++++++++ ++++++++   
 *
 * Then, there is an xover w/in locus 1, giving:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++---- -------- --------   
 *
 * Then, 1 x-over event b/w loci 1 and 2:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++---- ++++++++ ++++++++   
 *
 * No recombination events w/in locus 2:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++---- ++++++++ ++++++++   
 *
 * Two recombination events b/w loci 2 and 3,
 * which means double x-over, which means 
 * no change in the system:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++---- ++++++++ ++++++++   
 *
 * Finally, recombination w/in locus 3:
 *
 * 0      1 1      2 2      3 3      4
 * ++++++++ ++++---- ++++++++ ++++----   
 *
 * Thus, the parent's first gamete gives 
 * the following interval's to the offspring:
 *
 * [0, 1.5)
 * [2, 3.5)
 *
 * The parent's second gamete passed on the following
 * intervals:
 *
 * [1.5, 2.0)
 * [3.5, 4.0)
 *
 * Thus, the breakpoints we need to send for tree 
 * sequence recording are:
 *
 * [1.5, 2.0, 3.5, DBL_MAX],
 *
 * assuming that 4 is set as the "genome length"
 * of a table_collection.
 *
 * These values are stored in "expected_breakpoints"
 * in the fixture.
 *
 *
 * For the purpose of generating gametes, that "2.0"
 * breakpoint never needs to be stored, because parent 
 * gamete swapping, etc., takes care of that implicitly.
 *
 */

