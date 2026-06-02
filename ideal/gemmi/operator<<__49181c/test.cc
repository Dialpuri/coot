#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <sstream>
#include <vector>
#include <string>
#include "function.hh"

TEST(OracleTest, parallel_planes_t_stream_operator) {
    // Case 1: Two real residues from the PDB (chain A, residues 10 and 15)
    {
        coot::residue_spec_t rs1("A", 10);
        coot::residue_spec_t rs2("A", 15);
        std::vector<std::string> empty_names;
        coot::parallel_planes_t pp(rs1, rs2, empty_names, empty_names, "", "");

        std::ostringstream os;
        coot::stream_parallel_planes_t_gemmi(os, pp);
        std::string result = os.str();

        EXPECT_EQ(result, "pp-restr: [spec: mmdb::MinInt4 \"A\" 10 \"\"] [spec: mmdb::MinInt4 \"A\" 15 \"\"]");
    }

    // Case 2: Different chains (A and B)
    {
        coot::residue_spec_t rs1("A", 1);
        coot::residue_spec_t rs2("B", 1);
        std::vector<std::string> empty_names;
        coot::parallel_planes_t pp(rs1, rs2, empty_names, empty_names, "", "");

        std::ostringstream os;
        coot::stream_parallel_planes_t_gemmi(os, pp);
        std::string result = os.str();

        EXPECT_EQ(result, "pp-restr: [spec: mmdb::MinInt4 \"A\" 1 \"\"] [spec: mmdb::MinInt4 \"B\" 1 \"\"]");
    }

    // Case 3: Default/empty residue specs (edge case)
    {
        coot::residue_spec_t rs1;
        coot::residue_spec_t rs2;
        std::vector<std::string> empty_names;
        coot::parallel_planes_t pp(rs1, rs2, empty_names, empty_names, "", "");

        std::ostringstream os;
        coot::stream_parallel_planes_t_gemmi(os, pp);
        std::string result = os.str();

        EXPECT_EQ(result, "pp-restr: {residue-spec-not-set} {residue-spec-not-set}");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}