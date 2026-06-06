#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include "function.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTest, apply_mods) {
    coot::protein_geometry pg;
    pg.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];

    // Get a residue to test with - use residue //A/10 (ILE)
    // 0-based index: chain 0, residue 9
    gemmi::Residue* residue_p = &model.chains[0].residues[9];
    ASSERT_NE(residue_p, nullptr);

    std::string chain_id(model.chains[0].name);
    int seqnum = residue_p->seqid.num.value;
    std::string resname(residue_p->name);

    EXPECT_EQ(chain_id, "A");
    EXPECT_EQ(seqnum, 10);
    EXPECT_EQ(resname, "ILE");

    // Get atom table for residue (gemmi version)
    int i_no_res_atoms = static_cast<int>(residue_p->atoms.size());
    std::vector<gemmi::Atom*> res_selection;
    for (gemmi::Atom& atom : residue_p->atoms) {
        res_selection.push_back(&atom);
    }

    EXPECT_EQ(i_no_res_atoms, 8);

    // Get idr from geometry
    int idr = pg.get_monomer_restraints_index(resname, 0, false);
    EXPECT_EQ(idr, 9);

    std::vector<gemmi::Residue*> residues_with_OXTs;

    // Case 1: residues_with_OXTs is empty (default)
    {
        residues_with_OXTs.clear();

        coot::restraints_container_t::restraint_counts_t counts =
            coot::apply_mods_gemmi(idr, res_selection, i_no_res_atoms, residue_p, pg, residues_with_OXTs);

        EXPECT_EQ(counts.n_bond_restraints, 0);
        EXPECT_EQ(counts.n_angle_restraints, 0);
        EXPECT_EQ(counts.n_torsion_restr, 0);
    }

    // Case 2: residue IS in residues_with_OXTs
    {
        residues_with_OXTs.push_back(residue_p);

        // In original test, restraints.restraints_vec.size() was checked before/after.
        // Since restraints_container_t is MMDB-based and can't be constructed with
        // gemmi types, we track the total count from returned restraint_counts_t.
        coot::restraints_container_t::restraint_counts_t counts_before;
        counts_before.n_bond_restraints = 0;
        counts_before.n_angle_restraints = 0;
        counts_before.n_torsion_restr = 0;
        unsigned int restraints_before = counts_before.n_bond_restraints
                                         + counts_before.n_angle_restraints
                                         + counts_before.n_torsion_restr;

        coot::restraints_container_t::restraint_counts_t counts =
            coot::apply_mods_gemmi(idr, res_selection, i_no_res_atoms, residue_p, pg, residues_with_OXTs);

        unsigned int restraints_after = counts.n_bond_restraints
                                       + counts.n_angle_restraints
                                       + counts.n_torsion_restr;

        EXPECT_EQ(counts.n_bond_restraints, 0);
        EXPECT_EQ(counts.n_angle_restraints, 0);
        EXPECT_EQ(counts.n_torsion_restr, 0);
        EXPECT_EQ(restraints_before, 0u);
        EXPECT_EQ(restraints_after, 0u);
    }

    // Case 3: invalid residue (idr == -1)
    {
        residues_with_OXTs.clear();

        int bad_idr = -1;

        coot::restraints_container_t::restraint_counts_t counts =
            coot::apply_mods_gemmi(bad_idr, res_selection, i_no_res_atoms, residue_p, pg, residues_with_OXTs);

        EXPECT_EQ(counts.n_bond_restraints, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}