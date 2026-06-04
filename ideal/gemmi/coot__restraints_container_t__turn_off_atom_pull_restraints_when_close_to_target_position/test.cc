#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, turn_off_atom_pull_restraints_when_close_to_target_position) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model& model = st.models[0];
    ASSERT_FALSE(model.chains.empty());
    gemmi::Chain& chain = model.chains[0];
    ASSERT_GE(chain.residues.size(), 3u);

    // Get first atom from first 3 residues
    const gemmi::Atom* atom1 = &chain.residues[0].atoms[0];
    const gemmi::Atom* atom2 = &chain.residues[1].atoms[0];
    const gemmi::Atom* atom3 = &chain.residues[2].atoms[0];

    // Build atoms vector (indexed by position)
    std::vector<const gemmi::Atom*> atoms = {atom1, atom2, atom3};

    // Helper to make atom_spec_t from gemmi atom data
    // constructor: atom_spec_t(chain, resno, insertion_code, atom_name, alt_conf)
    auto make_spec = [&](const gemmi::Atom* at, const gemmi::Residue* res, const gemmi::Chain* ch) -> coot::atom_spec_t {
        std::string icode(1, res->seqid.icode);
        if (icode == " ") icode = "";
        return coot::atom_spec_t(ch->name, res->seqid.num.value, icode, at->name, "");
    };

    coot::atom_spec_t spec1 = make_spec(atom1, &chain.residues[0], &chain);
    coot::atom_spec_t spec2 = make_spec(atom2, &chain.residues[1], &chain);
    coot::atom_spec_t spec3 = make_spec(atom3, &chain.residues[2], &chain);

    // Use std::vector<simple_restraint> directly (free function version)
    std::vector<coot::simple_restraint> restraints;

    clipper::Coord_orth atom1_pos(atom1->pos.x, atom1->pos.y, atom1->pos.z);
    clipper::Coord_orth atom2_pos(atom2->pos.x, atom2->pos.y, atom2->pos.z);
    clipper::Coord_orth atom3_pos(atom3->pos.x, atom3->pos.y, atom3->pos.z);

    // Case 1: Add pull restraint with target VERY CLOSE to atom1 position (should be turned off)
    {
        coot::simple_restraint sr;
        sr.restraint_type = coot::restraint_type_t(256); // TARGET_POS_RESTRAINT
        sr.atom_spec = spec1;
        sr.atom_pull_target_pos = atom1_pos;
        sr.atom_index_1 = 0;
        sr.is_closed = false;
        restraints.push_back(sr);

        coot::atom_spec_t dragged_spec = spec2;

        std::vector<coot::atom_spec_t> result =
            coot::turn_off_atom_pull_restraints_when_close_to_target_position_gemmi(restraints, atoms, dragged_spec);

        EXPECT_EQ(result.size(), 1u);
        EXPECT_EQ(result[0], spec1);
    }

    // Case 2: Add pull restraint with target FAR from atom2 (should NOT be turned off)
    {
        clipper::Coord_orth far_target = atom2_pos + clipper::Coord_orth(10.0, 0.0, 0.0);
        coot::simple_restraint sr;
        sr.restraint_type = coot::restraint_type_t(256); // TARGET_POS_RESTRAINT
        sr.atom_spec = spec2;
        sr.atom_pull_target_pos = far_target;
        sr.atom_index_1 = 1;
        sr.is_closed = false;
        restraints.push_back(sr);

        coot::atom_spec_t dragged_spec = spec1;

        std::vector<coot::atom_spec_t> result =
            coot::turn_off_atom_pull_restraints_when_close_to_target_position_gemmi(restraints, atoms, dragged_spec);

        EXPECT_EQ(result.size(), 0u);
    }

    // Case 3: Add pull restraint with target close, but atom IS the dragged atom (should NOT be turned off)
    {
        coot::simple_restraint sr;
        sr.restraint_type = coot::restraint_type_t(256); // TARGET_POS_RESTRAINT
        sr.atom_spec = spec3;
        sr.atom_pull_target_pos = atom3_pos;
        sr.atom_index_1 = 2;
        sr.is_closed = false;
        restraints.push_back(sr);

        coot::atom_spec_t dragged_spec = spec3;

        std::vector<coot::atom_spec_t> result =
            coot::turn_off_atom_pull_restraints_when_close_to_target_position_gemmi(restraints, atoms, dragged_spec);

        EXPECT_EQ(result.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}