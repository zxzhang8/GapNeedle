#include "gapneedle/fasta_io.hpp"
#include "gapneedle/mapping_service.hpp"
#include "gapneedle/paf.hpp"

#include <cassert>
#include <fstream>
#include <iostream>

int main() {
  {
    const std::string seq = "ACGTN";
    const std::string rc = gapneedle::reverseComplement(seq);
    assert(rc == "NACGT");
  }

  {
    std::ofstream paf("/tmp/gapneedle_test.paf");
    paf << "q1\t100\t10\t40\t+\tt1\t120\t20\t50\t25\t30\t60\tcg:Z:30M\n";
    paf.close();

    auto recs = gapneedle::parsePaf("/tmp/gapneedle_test.paf", "t1", "q1");
    assert(recs.size() == 1);
    auto m = gapneedle::mapQueryToTargetDetail(recs[0], 15);
    assert(m.reason == "ok");
    assert(m.tPos.has_value());
    assert(m.tPos.value() == 25);
  }

  std::cout << "All tests passed\n";
  return 0;
}
