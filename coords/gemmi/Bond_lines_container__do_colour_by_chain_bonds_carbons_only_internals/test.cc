#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::string trim(const std::string &s) {
    size_t start = s.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

static int count_bonds(const Bond_lines_container_for_carbons &bonds) {
    int total = 0;
    for (const auto &bl : bonds.bonds) {
        total += static_cast<int>(bl.size());
    }
    return total;
}

TEST(OracleTest, do_colour_by_chain_bonds_carbons_only_internals_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];

    // residue 1 (MET)
    gemmi::Residue& res1 = chain.residues[0];
    // residue 2 (LYS)
    gemmi::Residue& res2 = chain.residues[1];

    // Find atoms in residue 1
    const gemmi::Atom* at_CA = nullptr;
    const gemmi::Atom* at_C  = nullptr;
    const gemmi::Atom* at_N  = nullptr;
    const gemmi::Atom* at_H  = nullptr;

    for (const gemmi::Atom& a : res1.atoms) {
        std::string aname = a.name;
        // Trim whitespace
        size_t start = aname.find_first_not_of(" \t");
        if (start == std::string::npos) continue;
        size_t end = aname.find_last_not_of(" \t");
        aname = aname.substr(start, end - start + 1);
        if (aname == "CA") at_CA = &a;
        if (aname == "C")  at_C  = &a;
        if (aname == "N")  at_N  = &a;
        if (aname == "H")  at_H  = &a;
    }

    // Find atoms in residue 2
    const gemmi::Atom* at_CA2 = nullptr;
    const gemmi::Atom* at_N2  = nullptr;
    for (const gemmi::Atom& a : res2.atoms) {
        std::string aname = a.name;
        size_t start = aname.find_first_not_of(" \t");
        if (start == std::string::npos) continue;
        size_t end = aname.find_last_not_of(" \t");
        aname = aname.substr(start, end - start + 1);
        if (aname == "CA") at_CA2 = &a;
        if (aname == "N")  at_N2  = &a;
    }

    ASSERT_NE(at_CA, nullptr);
    ASSERT_NE(at_C, nullptr);
    ASSERT_NE(at_N, nullptr);
    ASSERT_NE(at_H, nullptr);
    ASSERT_NE(at_CA2, nullptr);
    ASSERT_NE(at_N2, nullptr);

    Bond_lines_container_for_carbons bonds;
    int chain_idx = 1;
    std::vector<std::pair<bool, gemmi::Residue*>> het_residues;
    int uddHnd = -1;
    int udd_user_defined_atom_colour_index_handle = -1;
    int atom_colour_type = 20; // coot::COLOUR_BY_CHAIN_C_ONLY

    // === Case 1: C-C bond (same element carbon) ===
    {
        int before = count_bonds(bonds);
        std::string element1 = at_CA->element.name();
        std::string element2 = at_C->element.name();
        int iat_1 = 0;
        int iat_2 = 1;

        coot::Cartesian atom_1(at_CA->pos.x, at_CA->pos.y, at_CA->pos.z);
        coot::Cartesian atom_2(at_C->pos.x, at_C->pos.y, at_C->pos.z);

        EXPECT_EQ(trim(at_CA->name), "CA");
        EXPECT_EQ(trim(at_C->name), "C");
        EXPECT_EQ(before, 0);

        do_colour_by_chain_bonds_carbons_only_internals_gemmi(
            bonds,
            0, 0, chain_idx,
            at_CA, at_C, iat_1, iat_2,
            &het_residues,
            element1, element2,
            atom_1, atom_2,
            atom_colour_type,
            uddHnd,
            udd_user_defined_atom_colour_index_handle
        );

        int after = count_bonds(bonds);
        EXPECT_EQ(after, 1);
    }

    // === Case 2: C-N bond (different elements, neither H) ===
    {
        int before = count_bonds(bonds);
        std::string element1 = at_CA2->element.name();
        std::string element2 = at_N2->element.name();

        coot::Cartesian a1(at_CA2->pos.x, at_CA2->pos.y, at_CA2->pos.z);
        coot::Cartesian a2(at_N2->pos.x, at_N2->pos.y, at_N2->pos.z);

        EXPECT_EQ(trim(at_CA2->name), "CA");
        EXPECT_EQ(trim(at_N2->name), "N");
        EXPECT_EQ(before, 1);

        do_colour_by_chain_bonds_carbons_only_internals_gemmi(
            bonds,
            0, 0, chain_idx,
            at_CA2, at_N2, 10, 11,
            &het_residues,
            element1, element2,
            a1, a2,
            atom_colour_type,
            uddHnd,
            udd_user_defined_atom_colour_index_handle
        );

        int after = count_bonds(bonds);
        EXPECT_EQ(after, 3);
    }

    // === Case 3: C-H bond (hydrogen) ===
    {
        int before = count_bonds(bonds);
        std::string element1 = at_CA->element.name();
        std::string element2 = at_H->element.name();

        coot::Cartesian a1(at_CA->pos.x, at_CA->pos.y, at_CA->pos.z);
        coot::Cartesian a2(at_H->pos.x, at_H->pos.y, at_H->pos.z);

        EXPECT_EQ(trim(at_CA->name), "CA");
        EXPECT_EQ(trim(at_H->name), "H");
        EXPECT_EQ(before, 3);

        do_colour_by_chain_bonds_carbons_only_internals_gemmi(
            bonds,
            0, 0, chain_idx,
            at_CA, at_H, 20, 21,
            &het_residues,
            element1, element2,
            a1, a2,
            atom_colour_type,
            uddHnd,
            udd_user_defined_atom_colour_index_handle
        );

        int after = count_bonds(bonds);
        EXPECT_EQ(after, 3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}