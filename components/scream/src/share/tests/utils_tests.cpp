#include <catch2/catch.hpp>

#include "share/util/scream_utils.hpp"
#include "share/util/scream_time_stamp.hpp"

TEST_CASE("field_layout") {
  using namespace scream;

  std::string A = "A";
  std::string B = "B";
  std::string C = "C";
  std::string D = "D";
  std::string E = "E";
  std::string F = "F";
  std::string G = "G";

  using LOLS_type = std::list<std::list<std::string>>;

  // These three lists do not allow a superset from which they can all be
  // contiguously subviewed.
  LOLS_type lol1 = { {A,B}, {B,C}, {A,C} };
  REQUIRE(contiguous_superset(lol1).size()==0);

  // Input inner lists are not sorted
  REQUIRE_THROWS(contiguous_superset(LOLS_type{ {B,A} }));

  // The following should both allow the superset (A,B,C,D,E,F,G)
  // Note: lol3 is simply a shuffled version of lol2
  LOLS_type lol2 = { {A,B,C}, {B,C,D,E}, {C,D}, {C,D,E,F}, {D,E,F,G} };
  LOLS_type lol3 = { {D,E,F,G}, {C,D,E,F}, {A,B,C}, {C,D}, {B,C,D,E} };

  // Flipping a list is still a valid solution, so consider both tgt and its reverse.
  std::list<std::string> tgt = {A,B,C,D,E,F,G};
  std::list<std::string> tgt_rev = tgt;
  tgt_rev.reverse();

  auto superset2 = contiguous_superset(lol2);
  auto superset3 = contiguous_superset(lol3);
  REQUIRE ( (superset2==tgt || superset2==tgt_rev) );
  REQUIRE ( (superset3==tgt || superset3==tgt_rev) );
}

TEST_CASE ("time_stamp") {
  using namespace scream;
  using TS = util::TimeStamp;

  TS ts1 (2021,10,12,17,8,30);
  REQUIRE (ts1.get_years()==2021);
  REQUIRE (ts1.get_months()==10);
  REQUIRE (ts1.get_days()==12);
  REQUIRE (ts1.get_hours()==17);
  REQUIRE (ts1.get_minutes()==8);
  REQUIRE (ts1.get_seconds()==30);

  // Julian day = day_of_year.fraction_of_day, with day_of_year=0 at Jan 1st.
  REQUIRE (ts1.get_julian_day()==(284 + (17*3600+8*60+30)/86400.0));

  // Comparisons
  TS ts2 = ts1;
  ts2 += 10;
  REQUIRE (ts1<ts2);
  REQUIRE (ts2<=ts2);
  REQUIRE (ts2==ts2);

  // Cannot rewind time
  REQUIRE_THROWS (ts2+=-10);

}
