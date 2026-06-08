#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "ideal/simple-restraint.hh"
#include "function.hh"

TEST(OracleTest, restraints_container_t_add_chirals) {
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];

    auto find_residue = [&](const std::string& chain_name, int seq_num)
        -> std::pair<gemmi::Chain*, gemmi::Residue*>
    {
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != chain_name) continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return {&chain, &res};
                }
            }
        }
        return {nullptr, nullptr};
    };

    // Case 1: SER at chain A, residue 46 — has chiral center
    {
        auto cr = find_residue("A", 46);
        gemmi::Chain* chain_p  = cr.first;
        gemmi::Residue* res_p  = cr.second;
        ASSERT_NE(res_p, nullptr);

        std::string pdb_resname(res_p->name);
        EXPECT_EQ(pdb_resname, "SER");
        EXPECT_EQ(chain_p->name, "A");
        EXPECT_EQ(res_p->seqid.num.value, 46);

        int idr = geom.get_monomer_restraints_index(pdb_resname, 0, false);
        EXPECT_EQ(idr, 16);

        std::vector<gemmi::Atom*> res_selection;
        for (gemmi::Atom& atom : res_p->atoms) {
            res_selection.push_back(&atom);
        }
        int i_no_res_atoms = static_cast<int>(res_selection.size());
        EXPECT_EQ(i_no_res_atoms, 13);

        gemmi::CRA cra{chain_p, res_p, nullptr};
        std::vector<coot::simple_restraint> restraints_vec;

        EXPECT_EQ(restraints_vec.size(), 0u);

        int n_chirals = coot::add_chirals_gemmi(
            idr, res_selection, i_no_res_atoms, cra, geom, restraints_vec);

        EXPECT_EQ(restraints_vec.size(), 1u);
        EXPECT_EQ(n_chirals, 1);
    }

    // Case 2: GLY at chain A, residue 11 — no chiral center
    {
        auto cr2 = find_residue("A", 11);
        gemmi::Chain* chain_p2  = cr2.first;
        gemmi::Residue* res_p2  = cr2.second;
        ASSERT_NE(res_p2, nullptr);

        std::string pdb_resname2(res_p2->name);
        EXPECT_EQ(pdb_resname2, "GLY");
        EXPECT_EQ(chain_p2->name, "A");
        EXPECT_EQ(res_p2->seqid.num.value, 11);

        int idr2 = geom.get_monomer_restraints_index(pdb_resname2, 0, false);
        EXPECT_EQ(idr2, 5);

        std::vector<gemmi::Atom*> res_selection2;
        for (gemmi::Atom& atom : res_p2->atoms) {
            res_selection2.push_back(&atom);
        }
        int i_no_res_atoms2 = static_cast<int>(res_selection2.size());
        EXPECT_EQ(i_no_res_atoms2, 7);

        gemmi::CRA cra2{chain_p2, res_p2, nullptr};
        std::vector<coot::simple_restraint> restraints_vec2;

        size_t before2 = restraints_vec2.size();
        EXPECT_EQ(before2, 0u);

        int n_chirals2 = coot::add_chirals_gemmi(
            idr2, res_selection2, i_no_res_atoms2, cra2, geom, restraints_vec2);

        size_t after2 = restraints_vec2.size();
        EXPECT_EQ(after2, 0u);
        EXPECT_EQ(n_chirals2, 0);
    }

    // Case 3: ILE at chain A, residue 49 — another chiral residue
    {
        auto cr3 = find_residue("A", 49);
        gemmi::Chain* chain_p3  = cr3.first;
        gemmi::Residue* res_p3  = cr3.second;
        ASSERT_NE(res_p3, nullptr);

        std::string pdb_resname3(res_p3->name);
        EXPECT_EQ(pdb_resname3, "ILE");
        EXPECT_EQ(chain_p3->name, "A");
        EXPECT_EQ(res_p3->seqid.num.value, 49);

        int idr3 = geom.get_monomer_restraints_index(pdb_resname3, 0, false);
        EXPECT_GE(idr3, 0);

        std::vector<gemmi::Atom*> res_selection3;
        for (gemmi::Atom& atom : res_p3->atoms) {
            res_selection3.push_back(&atom);
        }
        int i_no_res_atoms3 = static_cast<int>(res_selection3.size());
        EXPECT_GT(i_no_res_atoms3, 0);

        gemmi::CRA cra3{chain_p3, res_p3, nullptr};
        std::vector<coot::simple_restraint> restraints_vec3;

        size_t before3 = restraints_vec3.size();
        EXPECT_EQ(before3, 0u);

        int n_chirals3 = coot::add_chirals_gemmi(
            idr3, res_selection3, i_no_res_atoms3, cra3, geom, restraints_vec3);

        EXPECT_GT(restraints_vec3.size(), before3);
        EXPECT_GT(n_chirals3, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}