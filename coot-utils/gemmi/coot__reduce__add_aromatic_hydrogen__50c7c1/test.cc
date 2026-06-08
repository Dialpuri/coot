#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <algorithm>
#include "function.hh"

// Helper: trim whitespace from both ends
static std::string trim(const std::string& s) {
    size_t first = s.find_first_not_of(" \t");
    if (first == std::string::npos) return "";
    size_t last  = s.find_last_not_of(" \t");
    return s.substr(first, last - first + 1);
}

// Helper: find a gemmi::Residue* by chain name and sequence number in the first model
static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int seqnum) {
    for (auto& chain : model.chains) {
        if (chain.name == chain_id) {
            for (auto& res : chain.residues) {
                if (res.seqid.num.value == seqnum) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

// Helper: delete an atom from a residue by its trimmed name
static bool delete_atom_from_residue(gemmi::Residue& res, const std::string& atom_name) {
    std::string trimmed = trim(atom_name);
    auto it = std::find_if(res.atoms.begin(), res.atoms.end(),
                           [&trimmed](const gemmi::Atom& a) { return trim(a.name) == trimmed; });
    if (it != res.atoms.end()) {
        res.atoms.erase(it);
        return true;
    }
    return false;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(AddAromaticHydrogenGemmiTest, valid_aromatic_h_addition) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];

    // Get PHE residue (seqnum 4, chain A)
    gemmi::Residue* phe_res = find_residue(model, "A", 4);
    ASSERT_NE(phe_res, nullptr);

    // Case 1: Valid aromatic H addition — remove HD1 first
    {
        delete_atom_from_residue(*phe_res, " HD1");
        int atoms_before = static_cast<int>(phe_res->atoms.size());

        EXPECT_EQ(atoms_before, 19);

        std::string H_at_name  = " HD1";
        std::string neighb_1   = " CG ";
        std::string neighb_2   = " CD1";
        std::string neighb_3   = " CE1";
        double bond_length     = 1.08;

        coot::reduce::add_aromatic_hydrogen_gemmi(H_at_name, neighb_1, neighb_2, neighb_3, bond_length, *phe_res);

        int atoms_after = static_cast<int>(phe_res->atoms.size());
        EXPECT_EQ(atoms_after, 20);
        EXPECT_EQ(atoms_after - atoms_before, 1);
    }

    // Case 2: Invalid neighbor — NONEXISTENT, should not add
    {
        int atoms_before = static_cast<int>(phe_res->atoms.size());

        EXPECT_EQ(atoms_before, 20);

        std::string H_at_name  = " HD1";
        std::string neighb_1   = " CG ";
        std::string neighb_2   = " NONEXISTENT";
        std::string neighb_3   = " CE1";
        double bond_length     = 1.08;

        coot::reduce::add_aromatic_hydrogen_gemmi(H_at_name, neighb_1, neighb_2, neighb_3, bond_length, *phe_res);

        int atoms_after = static_cast<int>(phe_res->atoms.size());
        EXPECT_EQ(atoms_after, 20);
        EXPECT_EQ(atoms_after - atoms_before, 0);
    }

    // Case 3: 4-arg overload with vector — remove HD2 first
    {
        delete_atom_from_residue(*phe_res, " HD2");
        int atoms_before = static_cast<int>(phe_res->atoms.size());

        EXPECT_EQ(atoms_before, 19);

        std::string H_at_name         = " HD2";
        std::string first_neigh       = " CG ";
        std::vector<std::string> second_neighb_vec = {" CD2", " CE2"};
        double bond_length            = 1.08;

        coot::reduce::add_aromatic_hydrogen_gemmi(H_at_name, first_neigh, second_neighb_vec, bond_length, phe_res);

        int atoms_after = static_cast<int>(phe_res->atoms.size());
        EXPECT_EQ(atoms_after, 20);
        EXPECT_EQ(atoms_after - atoms_before, 1);
    }
}
