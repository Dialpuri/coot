#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <tuple>
#include <string>
#include <iostream>
#include "function.hh"

// Helper to find a gemmi::CRA given chain_id, res_seq, icode, and atom_name
gemmi::CRA find_cra(const gemmi::Model& model,
                     const std::string& chain_id,
                     int res_seq,
                     char icode,
                     const std::string& atom_name)
{
    for (const gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (const gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == res_seq && res.seqid.icode == icode) {
                    for (const gemmi::Atom& atom : res.atoms) {
                        if (atom.name == atom_name) {
                            return {
                                const_cast<gemmi::Chain*>(&chain),
                                const_cast<gemmi::Residue*>(&res),
                                const_cast<gemmi::Atom*>(&atom)
                            };
                        }
                    }
                }
            }
        }
    }
    return {nullptr, nullptr, nullptr};
}

TEST(OracleTest, handle_long_bonded_atom) {
    // Load structure with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    EXPECT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];

    // Build Bond_lines_container_gemmi
    Bond_lines_container_gemmi blc;

    // Setup colour map — register chain "A"
    coot::my_atom_colour_map_t atom_colour_map;
    atom_colour_map.index_for_chain("A");

    // UDD handles — the gemmi function ignores these, pass original expected values
    int udd_handle_bond = 16777218;
    int udd_handle_atom_index = 16777217;
    int udd_user_defined_atom_colour_index_handle = 0;
    int atom_colour_type = 0;  // COLOUR_BY_CHAIN = 0

    // Model number is 0 (gemmi is 0-indexed)
    int model_number = 0;

    // Case 1: look up CA atom at //A/1/CA (the original MMDB test's fallback target)
    {
        gemmi::CRA cra = find_cra(model, "A", 1, ' ', "CA");
        EXPECT_NE(cra.atom, nullptr);

        // Verify atom name and element (gemmi stores unpadded names and elements)
        EXPECT_EQ(std::string(cra.atom->name), "CA");
        EXPECT_EQ(std::string(cra.atom->element.name()), "C");

        int total_before = 0;
        for (const auto& bl : blc.bonds) total_before += bl.size();

        // The original expected 5384 bonds before, 5392 after (delta 8).
        // Since we start with a fresh container (0 bonds), the delta is what matters.
        EXPECT_EQ(total_before, 0);

        int bonds_added = blc.handle_long_bonded_atom_gemmi(
            cra, model_number, udd_handle_bond, udd_handle_atom_index,
            udd_user_defined_atom_colour_index_handle,
            atom_colour_type, &atom_colour_map);

        int total_after = 0;
        for (const auto& bl : blc.bonds) total_after += bl.size();

        EXPECT_EQ(total_after - total_before, 8);
        EXPECT_EQ(bonds_added, 8);
    }

    // Case 2: look up CA atom at //A/10/CA (the original MMDB test's second fallback target)
    {
        gemmi::CRA cra2 = find_cra(model, "A", 10, ' ', "CA");
        EXPECT_NE(cra2.atom, nullptr);

        // Verify atom name and element
        EXPECT_EQ(std::string(cra2.atom->name), "CA");
        EXPECT_EQ(std::string(cra2.atom->element.name()), "C");

        int total_before2 = 0;
        for (const auto& bl : blc.bonds) total_before2 += bl.size();
        EXPECT_EQ(total_before2, 8);  // 8 bonds from case 1

        int bonds_added2 = blc.handle_long_bonded_atom_gemmi(
            cra2, model_number, udd_handle_bond, udd_handle_atom_index,
            udd_user_defined_atom_colour_index_handle,
            atom_colour_type, &atom_colour_map);

        int total_after2 = 0;
        for (const auto& bl : blc.bonds) total_after2 += bl.size();

        EXPECT_EQ(total_after2 - total_before2, 8);
        EXPECT_EQ(bonds_added2, 8);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
