#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/cif.hpp>
#include "function.hh"

TEST(OracleTest, pdbe_chem_comp_atom_depiction) {
    coot::protein_geometry geom;
    geom.init_standard();

    // Record state before
    std::size_t dict_size_before = geom.size();
    EXPECT_EQ(dict_size_before, std::size_t(41));

    // Build a cif::Loop with _pdbe_chem_comp_atom_depiction columns
    // Columns: comp_id, atom_id, element, model_Cartn_x, model_Cartn_y, pdbx_ordinal
    gemmi::cif::Loop loop;
    loop.tags.push_back("comp_id");
    loop.tags.push_back("atom_id");
    loop.tags.push_back("element");
    loop.tags.push_back("model_Cartn_x");
    loop.tags.push_back("model_Cartn_y");
    loop.tags.push_back("pdbx_ordinal");

    // Add row 1: ALA, CB, C, 1.0, 2.0, 1
    loop.add_row({"ALA", "CB", "C", "1.0", "2.0", "1"});
    // Add row 2: ALA, CG, C, 3.0, 4.0, 2
    loop.add_row({"ALA", "CG", "C", "3.0", "4.0", "2"});

    EXPECT_EQ(loop.length(), std::size_t(2));

    // IMOL_ENC_ANY = -2
    int imol_enc = -2;

    // Call the function
    coot::pdbe_chem_comp_atom_depiction_gemmi(geom, loop, imol_enc);

    // Record state after
    std::size_t dict_size_after = geom.size();
    EXPECT_EQ(dict_size_after, std::size_t(41));

    // Check if any depiction was added for ALA
    bool found_depiction = false;
    int atoms_in_depiction = 0;
    {
        int alaidx = geom.get_monomer_restraints_index("ALA", imol_enc, true);
        if (alaidx >= 0) {
            const auto &dep = geom[alaidx].second.depiction;
            if (!dep.empty()) {
                found_depiction = true;
                atoms_in_depiction = static_cast<int>(dep.atoms.size());
            }
        }
    }
    EXPECT_TRUE(found_depiction);
    EXPECT_EQ(atoms_in_depiction, 2);

    // CASE 2: Multiple comp_ids (should be rejected — guard requires comp_id_set.size() == 1)
    {
        gemmi::cif::Loop loop2;
        loop2.tags.push_back("comp_id");
        loop2.tags.push_back("atom_id");
        loop2.tags.push_back("element");
        loop2.tags.push_back("model_Cartn_x");
        loop2.tags.push_back("model_Cartn_y");
        loop2.tags.push_back("pdbx_ordinal");

        loop2.add_row({"ALA", "CB", "C", "1.0", "2.0", "1"});
        loop2.add_row({"GLY", "CA", "C", "3.0", "4.0", "2"});

        EXPECT_EQ(loop2.length(), std::size_t(2));

        int dep_before = atoms_in_depiction;
        coot::pdbe_chem_comp_atom_depiction_gemmi(geom, loop2, imol_enc);

        // Count ALA depiction atoms again
        int dep_after = 0;
        {
            int alaidx = geom.get_monomer_restraints_index("ALA", imol_enc, true);
            if (alaidx >= 0) {
                const auto &dep = geom[alaidx].second.depiction;
                dep_after = static_cast<int>(dep.atoms.size());
            }
        }
        EXPECT_EQ(dep_before, dep_after);
    }

    // CASE 3: Empty loop (no rows — should be a no-op)
    {
        gemmi::cif::Loop loop3;
        loop3.tags.push_back("comp_id");
        loop3.tags.push_back("atom_id");
        loop3.tags.push_back("element");
        loop3.tags.push_back("model_Cartn_x");
        loop3.tags.push_back("model_Cartn_y");
        loop3.tags.push_back("pdbx_ordinal");

        EXPECT_EQ(loop3.length(), std::size_t(0));

        EXPECT_NO_THROW(
            coot::pdbe_chem_comp_atom_depiction_gemmi(geom, loop3, imol_enc)
        );
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}