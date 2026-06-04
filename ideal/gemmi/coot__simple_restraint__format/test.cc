#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include "function.hh"

using restraint_type_t = coot::restraint_type_t;

TEST(simple_restraint_format, bond_positive) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];

    std::vector<coot::atom_spec_t> atom_vec;
    for (int r = 0; r < 3; r++) {
        if (static_cast<std::size_t>(r) >= chain.residues.size()) break;
        gemmi::Residue& res = chain.residues[r];
        for (std::size_t a = 0; a < res.atoms.size(); a++) {
            gemmi::Atom& at = res.atoms[a];
            // Build atom_spec_t from gemmi data using the 5-arg constructor
            coot::atom_spec_t spec(
                chain.name,                                      // chain_id
                res.seqid.num.value,                              // resno
                (res.seqid.icode == ' ' ? "" : std::string(1, res.seqid.icode)), // ins_code
                at.name,                                          // atom_name (padded)
                (at.altloc == '\0' ? "" : std::string(1, at.altloc))            // alt_conf
            );
            atom_vec.push_back(spec);
        }
    }

    ASSERT_GE(atom_vec.size(), 4u);

    // Case 1: BOND_RESTRAINT between atom 0 and atom 1
    {
        coot::simple_restraint r;
        r.restraint_type = restraint_type_t::BOND_RESTRAINT;
        r.atom_index_1 = 0;
        r.atom_index_2 = 1;
        std::string result = coot::format_gemmi(r, atom_vec, 4.0);
        EXPECT_EQ(result, "Bond N/1/A CA/1/A   z = 2.00");
    }

    // Case 2: BOND_RESTRAINT with negative distortion
    {
        coot::simple_restraint r;
        r.restraint_type = restraint_type_t::BOND_RESTRAINT;
        r.atom_index_1 = 0;
        r.atom_index_2 = 1;
        std::string result = coot::format_gemmi(r, atom_vec, -1.0);
        EXPECT_EQ(result, "Bond N/1/A CA/1/A ");
    }

    // Case 3: ANGLE_RESTRAINT
    {
        coot::simple_restraint r;
        r.restraint_type = restraint_type_t::ANGLE_RESTRAINT;
        r.atom_index_1 = 0;
        r.atom_index_2 = 1;
        r.atom_index_3 = 2;
        std::string result = coot::format_gemmi(r, atom_vec, 9.0);
        EXPECT_EQ(result, "Angle N/1/A CA/1/A C/1/A   z = 3.00");
    }

    // Case 4: TORSION_RESTRAINT
    {
        coot::simple_restraint r;
        r.restraint_type = restraint_type_t::TORSION_RESTRAINT;
        r.atom_index_1 = 0;
        r.atom_index_2 = 1;
        r.atom_index_3 = 2;
        r.atom_index_4 = 3;
        std::string result = coot::format_gemmi(r, atom_vec, 16.0);
        EXPECT_EQ(result, "Torsion N/1/A CA/1/A C/1/A O/1/A   z = 4.00");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}