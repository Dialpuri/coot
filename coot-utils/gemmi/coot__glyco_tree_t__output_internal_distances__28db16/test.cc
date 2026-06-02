#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

// Helper to read all lines from file
static std::vector<std::string> read_lines(const std::string &filename) {
    std::vector<std::string> lines;
    std::ifstream fin(filename.c_str());
    std::string line;
    while (std::getline(fin, line)) {
        lines.push_back(line);
    }
    return lines;
}

TEST(OracleTest, output_internal_distances) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    EXPECT_GE(st.models.size(), 1u);

    gemmi::Model &model = st.models[0];
    EXPECT_GE(model.chains.size(), 1u);

    gemmi::Chain &chain = model.chains[0];
    EXPECT_GE(chain.residues.size(), 1u);

    // Get first residue
    gemmi::Residue &residue = chain.residues[0];
    const std::string chain_id = chain.name;
    int seqnum = residue.seqid.num.value;
    const std::string resname = residue.name;

    // Count atoms in this residue
    size_t n_atoms = residue.atoms.size();

    EXPECT_EQ(chain_id, "A");
    EXPECT_EQ(seqnum, 1);
    EXPECT_EQ(resname, "MET");
    EXPECT_EQ(n_atoms, 19u);

    // We don't need a glyco_tree_t instance since the gemmi port is a free function
    // that doesn't use any member state

    double dist_crit = 5.0;

    // Case 1: self distances only (parent_p = nullptr)
    {
        std::string tmpfile = "/tmp/oracle_dist_case1.txt";
        std::ofstream f(tmpfile.c_str());
        if (f) {
            coot::output_internal_distances_gemmi(
                residue, nullptr, dist_crit, f, 1, chain_id);
            f.close();

            auto lines = read_lines(tmpfile);
            EXPECT_EQ(lines.size(), 27u);

            // Check first line
            EXPECT_NE(lines[0].find("[spec: model 1 \"A\"    1 \"\" \" N  \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[0].find("[spec: model 1 \"A\"    1 \"\" \" CA \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[0].find("1.48864"), std::string::npos);

            // Check last line (27th line, index 26)
            EXPECT_NE(lines[26].find("[spec: model 1 \"A\"    1 \"\" \" SD \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[26].find("[spec: model 1 \"A\"    1 \"\" \" CE \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[26].find("1.79326"), std::string::npos);
        }
    }

    // Case 2: self + parent distances (with a real parent residue)
    if (chain.residues.size() >= 2) {
        gemmi::Residue &parent_res = chain.residues[1];
        const std::string parent_resname = parent_res.name;
        int parent_seqnum = parent_res.seqid.num.value;

        EXPECT_EQ(parent_resname, "GLU");
        EXPECT_EQ(parent_seqnum, 2);

        std::string tmpfile = "/tmp/oracle_dist_case2.txt";
        std::ofstream f(tmpfile.c_str());
        if (f) {
            coot::output_internal_distances_gemmi(
                residue, &parent_res, dist_crit, f, 1, chain_id);
            f.close();

            auto lines = read_lines(tmpfile);
            EXPECT_EQ(lines.size(), 54u);

            // Check first line
            EXPECT_NE(lines[0].find("[spec: model 1 \"A\"    1 \"\" \" N  \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[0].find("1.48864"), std::string::npos);

            // Check last line (54th line, index 53)
            EXPECT_NE(lines[53].find("[spec: model 1 \"A\"    1 \"\" \" CE \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[53].find("[spec: model 1 \"A\"    2 \"\" \" CA \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[53].find("4.49328"), std::string::npos);

            // Check that cross-residue distances are present (line 28, index 27 is first cross-residue distance)
            EXPECT_NE(lines[27].find("[spec: model 1 \"A\"    1 \"\" \" N  \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[27].find("[spec: model 1 \"A\"    2 \"\" \" N  \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[27].find("2.73041"), std::string::npos);
        }
    }

    // Case 3: very tight distance cutoff (should produce fewer or no distances)
    {
        std::string tmpfile = "/tmp/oracle_dist_case3.txt";
        std::ofstream f(tmpfile.c_str());
        if (f) {
            double tight_crit = 2.5;
            coot::output_internal_distances_gemmi(
                residue, nullptr, tight_crit, f, 1, chain_id);
            f.close();

            auto lines = read_lines(tmpfile);
            EXPECT_EQ(lines.size(), 10u);

            // Check first line
            EXPECT_NE(lines[0].find("[spec: model 1 \"A\"    1 \"\" \" N  \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[0].find("1.48864"), std::string::npos);

            // Check last line (10th line, index 9)
            EXPECT_NE(lines[9].find("[spec: model 1 \"A\"    1 \"\" \" SD \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[9].find("[spec: model 1 \"A\"    1 \"\" \" CE \" \"\"]"), std::string::npos);
            EXPECT_NE(lines[9].find("1.79326"), std::string::npos);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}