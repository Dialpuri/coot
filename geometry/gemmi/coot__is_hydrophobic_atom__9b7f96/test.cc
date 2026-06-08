#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static gemmi::CRA find_atom_by_cid(const std::string &cid, gemmi::Structure &st) {
    gemmi::CRA cra{nullptr, nullptr, nullptr};

    if (cid.empty())
        return cra;

    // Parse CID: "//chain/residue_num/atom_name"
    std::string rest = cid;
    if (rest.size() >= 2 && rest[0] == '/' && rest[1] == '/')
        rest = rest.substr(2);

    size_t pos1 = rest.find('/');
    if (pos1 == std::string::npos)
        return cra;

    std::string chain_id = rest.substr(0, pos1);
    std::string rest2 = rest.substr(pos1 + 1);

    size_t pos2 = rest2.find('/');
    std::string res_seq_num_str;
    std::string atom_name;

    if (pos2 == std::string::npos) {
        res_seq_num_str = rest2;
        atom_name = "";
    } else {
        res_seq_num_str = rest2.substr(0, pos2);
        atom_name = rest2.substr(pos2 + 1);
    }

    int res_seq_num = 0;
    try {
        res_seq_num = std::stoi(res_seq_num_str);
    } catch (...) {
        return cra;
    }

    for (auto &model : st.models) {
        for (auto &chain : model.chains) {
            if (chain.name != chain_id)
                continue;
            for (auto &residue : chain.residues) {
                if (residue.seqid.num.value != res_seq_num)
                    continue;
                cra.chain   = &chain;
                cra.residue = &residue;
                if (atom_name.empty()) {
                    if (!residue.atoms.empty())
                        cra.atom = &residue.atoms[0];
                } else {
                    for (auto &atom : residue.atoms) {
                        std::string trimmed = atom.name;
                        auto start = trimmed.find_first_not_of(' ');
                        if (start != std::string::npos)
                            trimmed = trimmed.substr(start);
                        while (!trimmed.empty() && trimmed.back() == ' ')
                            trimmed.pop_back();
                        if (trimmed == atom_name) {
                            cra.atom = &atom;
                            return cra;
                        }
                    }
                }
            }
        }
    }
    return cra;
}

TEST(OracleTest, is_hydrophobic_atom) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: PHE 4 CB — hydrophobic sidechain atom
    {
        std::string cid = "//A/4/CB";
        gemmi::CRA cra = find_atom_by_cid(cid, st);
        ASSERT_NE(cra.atom, nullptr);
        bool result = coot::is_hydrophobic_atom_gemmi(cra);
        EXPECT_TRUE(result);
    }

    // Case 2: PHE 4 CZ — hydrophobic aromatic ring atom
    {
        std::string cid = "//A/4/CZ";
        gemmi::CRA cra = find_atom_by_cid(cid, st);
        ASSERT_NE(cra.atom, nullptr);
        bool result = coot::is_hydrophobic_atom_gemmi(cra);
        EXPECT_TRUE(result);
    }

    // Case 3: GLY 11 CA — main chain atom of glycine
    {
        std::string cid = "//A/11/CA";
        gemmi::CRA cra = find_atom_by_cid(cid, st);
        ASSERT_NE(cra.atom, nullptr);
        bool result = coot::is_hydrophobic_atom_gemmi(cra);
        EXPECT_TRUE(result);
    }

    // Case 4: GLU 2 OE1 — charged/polar sidechain atom
    {
        std::string cid = "//A/2/OE1";
        gemmi::CRA cra = find_atom_by_cid(cid, st);
        ASSERT_NE(cra.atom, nullptr);
        bool result = coot::is_hydrophobic_atom_gemmi(cra);
        EXPECT_FALSE(result);
    }

    // Case 5: VAL 7 CG1 — hydrophobic aliphatic sidechain
    {
        std::string cid = "//A/7/CG1";
        gemmi::CRA cra = find_atom_by_cid(cid, st);
        ASSERT_NE(cra.atom, nullptr);
        bool result = coot::is_hydrophobic_atom_gemmi(cra);
        EXPECT_TRUE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
