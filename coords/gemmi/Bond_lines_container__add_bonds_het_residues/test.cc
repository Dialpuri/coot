#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "api/molecules-container.hh"
#include "coot-utils/atom-selection-container.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTest, add_bonds_het_residues) {
    // case: main case with het residues
    {
        molecules_container_t mc;
        int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
        ASSERT_GT(imol, -1) << "PDB load failed";

        mmdb::Manager *mol = mc.get_mol(imol);
        ASSERT_NE(mol, nullptr);

        // Create atom selection container
        atom_selection_container_t asc = make_asc(mol);

        // Load structure with gemmi to get gemmi::Residue pointers
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

        // Find het residues (non-standard amino acids and nucleic acids)
        // The original test uses !res->isAminoacid() && !res->isNucleotide()
        std::vector<std::pair<bool, gemmi::Residue *>> het_residues;

        if (!st.models.empty()) {
            gemmi::Model &model = st.models[0];
            for (gemmi::Chain &chain : model.chains) {
                for (gemmi::Residue &res : chain.residues) {
                    if (!is_std_aa_or_nuc(res.name)) {
                        het_residues.push_back(std::make_pair(true, &res));
                    }
                }
            }
        }

        EXPECT_EQ(het_residues.size(), 1u);

        // Setup geometry
        coot::protein_geometry geom;
        geom.init_standard();

        // Construct Bond_lines_container_gemmi (derived class for bonds storage)
        Bond_lines_container_gemmi bonds;

        // UDD handles — in the gemmi port, these are not really used
        int udd_atom_index_handle = 0;
        int udd_user_defined_atom_colour_index_handle = 0;

        EXPECT_EQ(udd_atom_index_handle, 0);
        EXPECT_EQ(udd_user_defined_atom_colour_index_handle, 0);

        // Count bonds before
        int bonds_before = 0;
        for (auto &bl : bonds.bonds) bonds_before += bl.size();
        EXPECT_EQ(bonds_before, 0);

        // Call add_bonds_het_residues (atom_colour_type=0 for COLOUR_BY_CHAIN)
        EXPECT_NO_THROW(
            add_bonds_het_residues_gemmi(bonds, het_residues, asc, imol, 0,
                                         0, -1,
                                         udd_atom_index_handle, udd_user_defined_atom_colour_index_handle,
                                         &geom, true, true);
        );

        // Count bonds after
        int bonds_after = 0;
        for (auto &bl : bonds.bonds) bonds_after += bl.size();
        EXPECT_EQ(bonds_after, 0);

        // case: empty het_residues (edge case - should be no-op)
        std::vector<std::pair<bool, gemmi::Residue *>> empty_hets;
        int bonds_before2 = bonds_after;

        EXPECT_EQ(empty_hets.size(), 0u);
        EXPECT_EQ(bonds_before2, 0);

        EXPECT_NO_THROW(
            add_bonds_het_residues_gemmi(bonds, empty_hets, asc, imol, 0,
                                         0, -1, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle,
                                         &geom, true, true);
        );

        int bonds_after2 = 0;
        for (auto &bl : bonds.bonds) bonds_after2 += bl.size();
        EXPECT_EQ(bonds_after2, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}