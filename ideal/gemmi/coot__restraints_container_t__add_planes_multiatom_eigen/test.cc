#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_planes_multiatom_eigen) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    coot::protein_geometry geom;
    geom.init_standard();

    // === CASE 1: Aromatic residue with plane restraints (PHE) ===
    {
        // Find chain A, residue with seqnum 4 (PHE)
        gemmi::Chain* chain_p = nullptr;
        gemmi::Residue* res_p = nullptr;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& res : chain.residues) {
                    if (chain.name == "A" && res.seqid.num.value == 4) {
                        chain_p = &chain;
                        res_p = &res;
                    }
                }
            }
        }

        EXPECT_NE(res_p, nullptr);
        std::string pdb_resname(res_p->name);
        int idr = geom.get_monomer_restraints_index(pdb_resname, 0, false);

        // Build CRA vector for all atoms in this residue
        std::vector<gemmi::CRA> cras;
        for (gemmi::Atom& atom : res_p->atoms) {
            cras.push_back(gemmi::CRA{chain_p, res_p, &atom});
        }

        int i_no_res_atoms = static_cast<int>(cras.size());

        EXPECT_EQ(chain_p->name, "A");
        EXPECT_EQ(4, res_p->seqid.num.value);
        EXPECT_EQ(7, idr);
        EXPECT_EQ(11, i_no_res_atoms);
        EXPECT_EQ(2U, geom[idr].second.plane_restraint.size());

        size_t n_before = 0;
        EXPECT_EQ(0U, n_before);

        int n_added = coot::add_planes_multiatom_eigen_gemmi(
            idr, cras, i_no_res_atoms, res_p, chain_p, geom);

        size_t n_after = n_before + static_cast<size_t>(n_added);
        EXPECT_EQ(1, n_added);
        EXPECT_EQ(1U, n_after);
    }

    // === CASE 2: Non-aromatic residue with no plane restraints (MET) ===
    {
        // Find chain A, residue with seqnum 1 (MET)
        gemmi::Chain* chain_p = nullptr;
        gemmi::Residue* res_p = nullptr;
        for (gemmi::Model& model : st.models) {
            for (gemmi::Chain& chain : model.chains) {
                for (gemmi::Residue& res : chain.residues) {
                    if (chain.name == "A" && res.seqid.num.value == 1) {
                        chain_p = &chain;
                        res_p = &res;
                    }
                }
            }
        }

        EXPECT_NE(res_p, nullptr);
        std::string pdb_resname(res_p->name);
        int idr = geom.get_monomer_restraints_index(pdb_resname, 0, false);

        // Build CRA vector for all atoms in this residue
        std::vector<gemmi::CRA> cras;
        for (gemmi::Atom& atom : res_p->atoms) {
            cras.push_back(gemmi::CRA{chain_p, res_p, &atom});
        }

        int i_no_res_atoms = static_cast<int>(cras.size());

        EXPECT_EQ(chain_p->name, "A");
        EXPECT_EQ(1, res_p->seqid.num.value);
        EXPECT_EQ(12, idr);
        EXPECT_EQ(8, i_no_res_atoms);
        EXPECT_EQ(1U, geom[idr].second.plane_restraint.size());

        size_t n_before = 0;
        EXPECT_EQ(0U, n_before);

        int n_added = coot::add_planes_multiatom_eigen_gemmi(
            idr, cras, i_no_res_atoms, res_p, chain_p, geom);

        size_t n_after = n_before + static_cast<size_t>(n_added);
        EXPECT_EQ(0, n_added);
        EXPECT_EQ(0U, n_after);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}