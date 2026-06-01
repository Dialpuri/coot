#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "geometry/hb-types.hh"
#include "function.hh"

TEST(OracleTest, mark_donors_and_acceptors) {
    // Load molecule
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1);

    gemmi::Model& model = st.models[0];

    int imol = 0;

    // Init protein geometry (REQUIRED)
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Select all atoms in a single residue //A/10
    {
        auto atoms = coot::h_bonds::select_atoms_by_path(model, "//A/10");

        EXPECT_EQ(static_cast<int>(atoms.size()), 8);

        auto results = coot::h_bonds::mark_donors_and_acceptors_gemmi(atoms, atoms, geom, imol);

        int udd_handle = coot::h_bonds::mark_donors_and_acceptors_ud_handle;
        EXPECT_EQ(udd_handle, 16777217);

        // Atom 0: N  hb_type: 1 (HB_DONOR)
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms[0].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 1);
        }
        // Atom 1: CA hb_type: 0 (HB_NEITHER)
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms[1].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);
        }
        // Atom 2: C  hb_type: 0 (HB_NEITHER)
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms[2].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);
        }
        // Atom 3: O  hb_type: 2 (HB_ACCEPTOR)
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms[3].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 2);
        }
        // Atom 4: CB hb_type: 0 (HB_NEITHER)
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms[4].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);
        }
        // Atom 5: CG1 hb_type: 0 (HB_NEITHER)
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms[5].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);
        }
        // Atom 6: CG2 hb_type: 0 (HB_NEITHER)
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms[6].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);
        }
        // Atom 7: CD1 hb_type: 0 (HB_NEITHER)
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms[7].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);
        }
    }

    // Case 2: Two different residues //A/10 + //A/20
    {
        auto atoms1 = coot::h_bonds::select_atoms_by_path(model, "//A/10");
        auto atoms2 = coot::h_bonds::select_atoms_by_path(model, "//A/20");

        EXPECT_EQ(static_cast<int>(atoms1.size()), 8);
        EXPECT_EQ(static_cast<int>(atoms2.size()), 5);

        auto results = coot::h_bonds::mark_donors_and_acceptors_gemmi(atoms1, atoms2, geom, imol);

        int udd_handle = coot::h_bonds::mark_donors_and_acceptors_ud_handle;
        EXPECT_EQ(udd_handle, 16777217);

        // sel1 atoms: N=1, CA=0, C=0, O=2, CB=0, CG1=0, CG2=0, CD1=0
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms1[0].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 1);  // N
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms1[1].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);  // CA
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms1[2].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);  // C
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms1[3].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 2);  // O
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms1[4].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);  // CB
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms1[5].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);  // CG1
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms1[6].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);  // CG2
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms1[7].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);  // CD1
        }

        // sel2 atoms: N=1, CA=0, C=0, O=2, CB=0
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms2[0].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 1);  // N
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms2[1].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);  // CA
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms2[2].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);  // C
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms2[3].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 2);  // O
        }
        {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms2[4].atom, udd_handle, hb_type);
            EXPECT_EQ(hb_type, 0);  // CB
        }
    }

    // Case 3: Select whole chain A
    {
        auto atoms = coot::h_bonds::select_atoms_by_path(model, "//A");

        EXPECT_EQ(static_cast<int>(atoms.size()), 2106);

        auto results = coot::h_bonds::mark_donors_and_acceptors_gemmi(atoms, atoms, geom, imol);

        int udd_handle = coot::h_bonds::mark_donors_and_acceptors_ud_handle;
        EXPECT_EQ(udd_handle, 16777217);

        // Count atoms by hb_type
        int n_assigned = 0;
        int n_unassigned = 0;
        int n_donor = 0;
        int n_acceptor = 0;
        int n_both = 0;
        for (size_t i = 0; i < atoms.size(); i++) {
            int hb_type = -999;
            coot::h_bonds::get_hb_type(results, atoms[i].atom, udd_handle, hb_type);
            if (hb_type != -1) n_assigned++;  // HB_UNASSIGNED = -1
            else n_unassigned++;
            if (hb_type == 1) n_donor++;      // HB_DONOR = 1
            if (hb_type == 2) n_acceptor++;   // HB_ACCEPTOR = 2
            if (hb_type == 3) n_both++;       // HB_BOTH = 3
        }
        EXPECT_EQ(n_assigned, 2106);
        EXPECT_EQ(n_unassigned, 0);
        EXPECT_EQ(n_donor, 360);
        EXPECT_EQ(n_acceptor, 345);
        EXPECT_EQ(n_both, 37);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}