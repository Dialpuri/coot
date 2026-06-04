#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>  // read_structure_file
#include <gemmi/calculate.hpp>
#include <iostream>
#include "api/molecules-container.hh"
#include "api/coot-molecule.hh"
#include "function.hh"

TEST(OracleTest, sfcalc_genmap) {
    {
        // case: main flow - replicate oracle setup
        molecules_container_t mc;

        // Load PDB/CIF
        int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
        EXPECT_GE(imol, 0);

        // Load MTZ for map coefficients (FWT/PHWT)
        int mtz_imol = mc.read_mtz(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
            "FWT", "PHWT", "", false, false);
        // Oracle printed: INPUT  mtz_imol: 1
        EXPECT_EQ(mtz_imol, 1);

        coot::molecule_t& mol = mc[imol];
        coot::molecule_t& mtz_mol = mc[mtz_imol];

        // Get mmdb::Manager* from PDB molecule
        mmdb::Manager* mol_p = mol.atom_sel.mol;
        EXPECT_NE(mol_p, nullptr);

        // Try to get F_obs from MTZ molecule
        // fill_fobs_sigfobs() throws because MTZ only has FWT/PHWT, not F/SIGF
        const clipper::HKL_data<clipper::data32::F_sigF>* fobs_p = nullptr;
        try {
            mtz_mol.fill_fobs_sigfobs();
            fobs_p = mtz_mol.get_original_fobs_sigfobs();
        } catch (...) {
            fobs_p = nullptr;
        }

        // Oracle behavior: fill_fobs_sigfobs() fails → fobs_p is null
        EXPECT_EQ(fobs_p, nullptr);

        // Retry with F and SIGF columns
        if (!fobs_p || fobs_p->num_obs() == 0) {
            int mtz2 = -1;
            try {
                mtz2 = mc.read_mtz(
                    "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
                    "F", "SIGF", "", false, false);
            } catch (...) {
                mtz2 = -1;
            }
            if (mtz2 >= 0) {
                coot::molecule_t& mtz2_mol = mc[mtz2];
                try {
                    mtz2_mol.fill_fobs_sigfobs();
                    fobs_p = mtz2_mol.get_original_fobs_sigfobs();
                } catch (...) {
                    fobs_p = nullptr;
                }
            }
        }

        // Retry with F_obs/sigma_F columns
        if (!fobs_p || fobs_p->num_obs() == 0) {
            int mtz3 = -1;
            try {
                mtz3 = mc.read_mtz(
                    "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
                    "F_obs", "sigma_F", "", false, false);
            } catch (...) {
                mtz3 = -1;
            }
            if (mtz3 >= 0) {
                coot::molecule_t& mtz3_mol = mc[mtz3];
                try {
                    mtz3_mol.fill_fobs_sigfobs();
                    fobs_p = mtz3_mol.get_original_fobs_sigfobs();
                } catch (...) {
                    fobs_p = nullptr;
                }
            }
        }

        // Oracle behavior: after all retries, fobs_p is still null
        // Oracle would print: "OUTPUT no_fobs_final: Could not get F_obs from any column combo"
        // and then return 1 (exit early, never calling sfcalc_genmap)
        EXPECT_EQ(fobs_p, nullptr);
    }
}

TEST(OracleTest, sfcalc_genmap_gemmi_structure_load) {
    // Verify that the gemmi version can at least load the structure
    gemmi::Structure st = gemmi::read_structure_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    EXPECT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];
    // Count atoms manually since gemmi::Model has no n_atoms()
    size_t atom_count = 0;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            atom_count += res.atoms.size();
        }
    }
    EXPECT_GT(atom_count, 0u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}