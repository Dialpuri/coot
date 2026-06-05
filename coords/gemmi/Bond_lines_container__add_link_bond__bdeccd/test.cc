#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "coot-utils/atom-selection-container.hh"
#include "coords/Bond_lines.hh"
#include "api/molecules-container.hh"

using namespace coot;

TEST(OracleTest, add_link_bond) {
    molecules_container_t mc;
    mc.geometry_init_standard();
    int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(imol, 0);

    // Get the molecule and its gemmi structure
    mmdb::Manager *mol = mc.get_mol(imol);
    ASSERT_NE(mol, nullptr);

    // Get gemmi structure
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u);
    gemmi::Model& model = st.models[0];

    atom_selection_container_t asc = make_asc(mol);

    // Use gemmi helper class that exposes protected/private members
    Bond_lines_container_gemmi bonds;

    // Count total bonds helper
    auto total_bonds = [&]() {
        int total = 0;
        for (const auto &bl : bonds.bonds) total += bl.size();
        return total;
    };

    // Case 1: valid link between A/1/N and B/1/CA (cross-chain, no pre-existing bond)
    // Oracle printed bonds_before=0, bonds_after=0
    {
        gemmi_link_info_t link;
        link.chain1 = "A";
        link.seq_num1 = 1;
        link.res_name1 = "MET";
        link.atom_name1 = "N";
        link.icode1 = ' ';

        link.chain2 = "B";
        link.seq_num2 = 1;
        link.res_name2 = "MET";
        link.atom_name2 = "CA";
        link.icode2 = ' ';

        int total_before = total_bonds();
        add_link_bond_gemmi(bonds, model, -1, -1, 0, link);
        int total_after = total_bonds();

        EXPECT_EQ(total_before, 0);
        EXPECT_EQ(total_after, 0);
    }

    // Case 2: link to non-existent residue number
    // Oracle printed bonds_before=0, bonds_after=0
    {
        gemmi_link_info_t link;
        link.chain1 = "A";
        link.seq_num1 = 9999;
        link.res_name1 = "MET";
        link.atom_name1 = "N";
        link.icode1 = ' ';

        link.chain2 = "A";
        link.seq_num2 = 9999;
        link.res_name2 = "MET";
        link.atom_name2 = "CA";
        link.icode2 = ' ';

        int total_before = total_bonds();
        add_link_bond_gemmi(bonds, model, -1, -1, 0, link);
        int total_after = total_bonds();

        EXPECT_EQ(total_before, 0);
        EXPECT_EQ(total_after, 0);
    }

    // Case 3 (complementary): Bond_lines_container default construction works
    // This verifies that the container is usable with default construction
    {
        Bond_lines_container_gemmi bonds_default;
        EXPECT_EQ(bonds_default.total_bond_count(), 0);
    }

    // Case 4 (complementary): add_link_bond with proper atom indices
    // Find actual atoms from the gemmi model for A/1/N and B/1/CA
    // Since N is nitrogen and CA is carbon, they have different elements,
    // so add_link_bond_templ creates TWO half-bonds (one per element colour)
    {
        // Find atoms in gemmi model
        gemmi::CRA cra1 = find_atom_cra_in_model(model, "A", 1, "MET", "N", ' ');
        gemmi::CRA cra2 = find_atom_cra_in_model(model, "B", 1, "MET", "CA", ' ');

        if (cra1.atom != nullptr && cra2.atom != nullptr) {
            gemmi_link_info_t link;
            link.chain1 = "A";
            link.seq_num1 = 1;
            link.res_name1 = "MET";
            link.atom_name1 = "N";
            link.icode1 = ' ';

            link.chain2 = "B";
            link.seq_num2 = 1;
            link.res_name2 = "MET";
            link.atom_name2 = "CA";
            link.icode2 = ' ';

            int total_before = total_bonds();
            add_link_bond_gemmi(bonds, model, -1, -1, 0, link);
            int total_after = total_bonds();
            EXPECT_GT(total_after, total_before);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}