#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <fstream>
#include <string>
#include <set>
#include <vector>
#include "function.hh"

// Helper to count total residues in a gemmi::Structure
static int count_residues(const gemmi::Structure& st) {
    int total = 0;
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            total += static_cast<int>(chain.residues.size());
        }
    }
    return total;
}

// Helper to count carbohydrate residues
static int count_carb_residues(const gemmi::Structure& st) {
    std::set<std::string> cho_set;
    cho_set.insert("NAG"); cho_set.insert("MAN"); cho_set.insert("BMA"); cho_set.insert("FUL");
    cho_set.insert("FUC"); cho_set.insert("XYP"); cho_set.insert("SIA"); cho_set.insert("GAL");
    cho_set.insert("NDG"); cho_set.insert("BGC"); cho_set.insert("A2G");

    int total = 0;
    for (const gemmi::Model& model : st.models) {
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& res : chain.residues) {
                if (cho_set.find(res.name) != cho_set.end()) {
                    total++;
                }
            }
        }
    }
    return total;
}

TEST(OracleTest, delete_all_carbohydrate) {
    // Case 1: Default PDB with no carbohydrates — should return false
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        int residues_before = count_residues(st);
        int carbs_before = count_carb_residues(st);

        EXPECT_EQ(residues_before, 268);
        EXPECT_EQ(carbs_before, 0);

        bool deleted = coot::util::delete_all_carbohydrate_gemmi(&st);

        int residues_after = count_residues(st);
        int carbs_after = count_carb_residues(st);

        EXPECT_FALSE(deleted);
        EXPECT_EQ(residues_after, 268);
        EXPECT_EQ(carbs_after, 0);
    }

    // Case 2: PDB with explicit carbohydrate residues (NAG, MAN)
    {
        std::string pdb_content =
            "HETATM    1  C1  NAG A  1       0.000   0.000   0.000  1.00  0.00           C  \n"
            "HETATM    2  C2  NAG A  1       1.000   0.000   0.000  1.00  0.00           C  \n"
            "HETATM    3  C3  NAG A  1       2.000   0.000   0.000  1.00  0.00           C  \n"
            "HETATM    4  O1  NAG A  1       3.000   0.000   0.000  1.00  0.00           O  \n"
            "HETATM    5  N   ALA A  2       4.000   0.000   0.000  1.00  0.00           N  \n"
            "HETATM    6  CA  ALA A  2       5.000   0.000   0.000  1.00  0.00           C  \n"
            "HETATM    7  C   ALA A  2       6.000   0.000   0.000  1.00  0.00           C  \n"
            "HETATM    8  O   ALA A  2       7.000   0.000   0.000  1.00  0.00           O  \n"
            "HETATM    9  C1  MAN A  3       8.000   0.000   0.000  1.00  0.00           C  \n"
            "END\n";

        std::string tmp_path = "/tmp/oracle_carb_test.pdb";
        std::ofstream ofs(tmp_path);
        ofs << pdb_content;
        ofs.close();

        gemmi::Structure st = gemmi::read_pdb_file(tmp_path);

        int residues_before = count_residues(st);
        int carbs_before = count_carb_residues(st);

        EXPECT_EQ(residues_before, 3);
        EXPECT_EQ(carbs_before, 2);

        bool deleted = coot::util::delete_all_carbohydrate_gemmi(&st);

        int residues_after = count_residues(st);
        int carbs_after = count_carb_residues(st);

        EXPECT_TRUE(deleted);
        EXPECT_EQ(residues_after, 1);
        EXPECT_EQ(carbs_after, 0);
    }

    // Case 3: nullptr input — should return false
    {
        bool deleted = coot::util::delete_all_carbohydrate_gemmi(nullptr);
        EXPECT_FALSE(deleted);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}