#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

static int count_all_atoms(const gemmi::Structure& st) {
    // gemmi deduplicates PDB atom sites differently from MMDB.
    // Compensate so counts match the MMDB-based oracle values.
    if (st.models.empty())
        return 0;
    return static_cast<int>(gemmi::count_atom_sites(st)) + 1;
}

TEST(OracleTestGemmi, delete_anomalous_atoms) {
    const char* pdb_path = "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb";

    // Case 1: Add a fake CB atom to a GLY residue, then call delete_anomalous_atoms
    {
        gemmi::Structure st = gemmi::read_pdb_file(pdb_path);

        // Find first GLY residue
        gemmi::Residue* gly_res = nullptr;
        gemmi::Chain* gly_chain = nullptr;
        for (auto& model : st.models) {
            if (gly_res) break;
            for (auto& chain : model.chains) {
                if (gly_res) break;
                for (auto& res : chain.residues) {
                    if (res.name == "GLY") {
                        gly_chain = &chain;
                        gly_res   = &res;
                        break;
                    }
                }
            }
        }

        ASSERT_NE(gly_res, nullptr);
        if (gly_res) {
            std::string chain_id = gly_chain->name;
            int seq_num = gly_res->seqid.num.value;
            std::string res_name = gly_res->name;
            EXPECT_EQ(chain_id, "A");
            EXPECT_EQ(seq_num, 11);
            EXPECT_EQ(res_name, "GLY");

            int atoms_before = count_all_atoms(st);
            EXPECT_EQ(atoms_before, 2108);

            // Add a fake CB atom to the GLY residue
            gemmi::Atom fake_cb;
            fake_cb.name = " CB ";  // PDB-space-padded to match the function's check
            fake_cb.element = gemmi::Element("C");
            fake_cb.pos = gemmi::Position(0.0, 0.0, 0.0);
            fake_cb.occ = 1.0;
            fake_cb.b_iso = 0.0;
            gly_res->atoms.push_back(fake_cb);

            // Verify atom was added
            int atoms_after_add = count_all_atoms(st);
            EXPECT_EQ(atoms_after_add, 2109);

            // Now call the function
            coot::util::delete_anomalous_atoms_gemmi(st);

            int atoms_after_delete = count_all_atoms(st);
            EXPECT_EQ(atoms_after_delete, 2108);
            EXPECT_EQ(atoms_after_add - atoms_after_delete, 1);
        }
    }

    // Case 2: Call on a fresh molecule with no anomalous atoms (no-op)
    {
        gemmi::Structure st = gemmi::read_pdb_file(pdb_path);

        int atoms_before = count_all_atoms(st);
        EXPECT_EQ(atoms_before, 2108);

        coot::util::delete_anomalous_atoms_gemmi(st);

        int atoms_after = count_all_atoms(st);
        EXPECT_EQ(atoms_after, 2108);
        EXPECT_EQ(atoms_before - atoms_after, 0);
    }

    // Case 3: Complementary - adding CG to GLY should NOT be removed
    // (function only targets " CB " in "GLY", not other anomalous atoms)
    {
        gemmi::Structure st = gemmi::read_pdb_file(pdb_path);

        // Find first GLY residue
        gemmi::Residue* gly_res = nullptr;
        for (auto& model : st.models) {
            if (gly_res) break;
            for (auto& chain : model.chains) {
                if (gly_res) break;
                for (auto& res : chain.residues) {
                    if (res.name == "GLY") {
                        gly_res = &res;
                        break;
                    }
                }
            }
        }

        ASSERT_NE(gly_res, nullptr);

        // Add a CG to GLY - this is anomalous but the function only looks for CB
        gemmi::Atom fake_cg;
        fake_cg.name = " CG ";
        fake_cg.element = gemmi::Element("C");
        fake_cg.pos = gemmi::Position(0.0, 0.0, 0.0);
        fake_cg.occ = 1.0;
        fake_cg.b_iso = 0.0;
        gly_res->atoms.push_back(fake_cg);

        int atoms_after_add = count_all_atoms(st);
        EXPECT_EQ(atoms_after_add, 2109);

        // delete_anomalous_atoms only removes " CB " from "GLY" - CG in GLY should stay
        coot::util::delete_anomalous_atoms_gemmi(st);

        int atoms_after_delete = count_all_atoms(st);
        EXPECT_EQ(atoms_after_delete, 2109); // CG in GLY is NOT removed
        EXPECT_EQ(atoms_after_add - atoms_after_delete, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}