#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>
#include "function.hh"

// Helper: find a gemmi::Atom by CID string like "//A/10/N"
static gemmi::Atom* find_atom_by_cid(const std::string &cid, gemmi::Structure &st) {
    if (cid.size() < 4 || cid[0] != '/' || cid[1] != '/')
        return nullptr;
    std::string rest = cid.substr(2);
    
    auto chain_end = rest.find('/');
    if (chain_end == std::string::npos) return nullptr;
    std::string chain_name = rest.substr(0, chain_end);
    
    auto res_end = rest.find('/', chain_end + 1);
    if (res_end == std::string::npos) return nullptr;
    std::string res_seq_str = rest.substr(chain_end + 1, res_end - chain_end - 1);
    
    std::string atom_name = rest.substr(res_end + 1);
    
    int res_seq = 0;
    try { res_seq = std::stoi(res_seq_str); } catch (...) { return nullptr; }
    
    for (gemmi::Model &model : st.models) {
        for (gemmi::Chain &chain : model.chains) {
            if (chain.name != chain_name) continue;
            for (gemmi::Residue &res : chain.residues) {
                if (res.seqid.num.value != res_seq) continue;
                for (gemmi::Atom &atom : res.atoms) {
                    if (atom.name == atom_name)
                        return &atom;
                }
            }
        }
    }
    return nullptr;
}

// Pad atom name to 4 characters to match MMDB GetAtomName() behavior
static std::string pad_atom_name(const std::string &name) {
    if (name.size() >= 4) return name.substr(0, 4);
    std::string padded = name;
    while (padded.size() < 4) padded += ' ';
    // 1-char and 2-char names are right-justified in PDB (leading space)
    if (name.size() <= 2 && name.size() >= 1) {
        padded.insert(0, 1, ' ');
        padded.pop_back();
    }
    return padded;
}

TEST(OracleTest, bonded_atom_pair_t) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: valid bonded_atom_pair with two real atoms (N and CA of residue A/10)
    {
        gemmi::Atom *atom_1 = find_atom_by_cid("//A/10/N", st);
        gemmi::Atom *atom_2 = find_atom_by_cid("//A/10/CA", st);
        double residual = 0.35;

        Bond_lines_container::Bond_lines_container::bonded_atom_pair_t pair =
            Bond_lines_container::Bond_lines_container::bonded_atom_pair_t::bonded_atom_pair_t_gemmi(atom_1, atom_2, residual);

        EXPECT_TRUE(pair.atom_1 != nullptr);
        EXPECT_TRUE(pair.atom_2 != nullptr);
        EXPECT_NEAR(pair.residual, 0.35, 1e-4);
        if (pair.atom_1) {
            EXPECT_EQ(pad_atom_name(pair.atom_1->name), " N  ");
        }
        if (pair.atom_2) {
            EXPECT_EQ(pad_atom_name(pair.atom_2->name), " CA ");
        }
    }

    // Case 2: bonded_atom_pair with one nullptr (missing atom)
    {
        gemmi::Atom *atom_1 = find_atom_by_cid("//A/9999/N", st);  // non-existent
        gemmi::Atom *atom_2 = find_atom_by_cid("//A/10/CA", st);
        double residual = 1.5;

        Bond_lines_container::Bond_lines_container::bonded_atom_pair_t pair =
            Bond_lines_container::Bond_lines_container::bonded_atom_pair_t::bonded_atom_pair_t_gemmi(atom_1, atom_2, residual);

        EXPECT_EQ(pair.atom_1, nullptr);
        EXPECT_TRUE(pair.atom_2 != nullptr);
        EXPECT_NEAR(pair.residual, 1.5, 1e-4);
    }

    // Case 3: bonded_atom_pair with both nullptrs
    {
        gemmi::Atom *atom_1 = find_atom_by_cid("//Z/999/N", st);
        gemmi::Atom *atom_2 = find_atom_by_cid("//Z/999/CA", st);
        double residual = 5.0;

        Bond_lines_container::Bond_lines_container::bonded_atom_pair_t pair =
            Bond_lines_container::Bond_lines_container::bonded_atom_pair_t::bonded_atom_pair_t_gemmi(atom_1, atom_2, residual);

        EXPECT_EQ(pair.atom_1, nullptr);
        EXPECT_EQ(pair.atom_2, nullptr);
        EXPECT_NEAR(pair.residual, 5.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}