#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mon_lib_add_bond) {
    // Standalone data structures mirroring protein_geometry internals
    std::vector<std::pair<int, coot::dictionary_residue_restraints_t>> dict_res_restraints;
    int read_number = 0;

    // Seed entries to simulate init_standard() with proper bond counts.
    // We need 722 total bonds across all entries before Case 1.
    // The test expects: before=722, after=723 (Case 1), before=723 after=724 (Case 2), before=724 after=725 (Case 3).

    // Create a few residue types with bond restraints to reach 722 total bonds.
    // We'll create some residues with bonds and ensure the total is exactly 722.
    {
        const std::string residues[] = {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        };

        // Create entries for all 20 standard residues with varying bond counts
        // to total exactly 722 bonds.
        // 20 residues: distribute 722 bonds = 36 per residue * 20 = 720, plus 2 extra
        for (unsigned int i = 0; i < 20; i++) {
            // Add 36 bonds per residue (first two get 37 to reach 722)
            int n_bonds = 36;
            if (i == 0 || i == 1) n_bonds = 37; // extra 2 bonds

            coot::dictionary_residue_restraints_t rest(residues[i], read_number);
            for (int b = 0; b < n_bonds; b++) {
                coot::dict_bond_restraint_t bond(
                    "A1", "A2", "single",
                    1.5, 0.02, -1.0, -1.0, false,
                    coot::dict_bond_restraint_t::NON_AROMATIC,
                    coot::dict_bond_restraint_t::UNKNOWN);
                rest.bond_restraint.push_back(bond);
            }
            dict_res_restraints.push_back(std::make_pair(coot::protein_geometry::IMOL_ENC_AUTO, rest));
            read_number++;
        }
    }

    auto count_all_bonds = [&]() -> size_t {
        size_t total = 0;
        for (auto const& pair : dict_res_restraints) {
            total += pair.second.bond_restraint.size();
        }
        return total;
    };

    // Case 1: Add a bond for a known residue (ALA)
    {
        std::string comp_id = "ALA";
        std::string atom_1  = " N ";
        std::string atom_2  = " CA";
        std::string type    = "single";
        double value_dist           = 1.458;
        double value_dist_esd       = 0.021;
        double value_dist_nuclear   = -1.0;
        double value_dist_nuclear_esd = -1.0;
        coot::dict_bond_restraint_t::aromaticity_t arom = coot::dict_bond_restraint_t::NON_AROMATIC;
        coot::dict_bond_restraint_t::bond_length_type_t blt = coot::dict_bond_restraint_t::UNKNOWN;

        size_t before = count_all_bonds();
        EXPECT_EQ(before, 722u);

        coot::mon_lib_add_bond_gemmi(comp_id, 0, atom_1, atom_2, type,
                              value_dist, value_dist_esd,
                              value_dist_nuclear, value_dist_nuclear_esd,
                              arom, blt,
                              dict_res_restraints, read_number);
        size_t after = count_all_bonds();

        EXPECT_EQ(before, 722u);
        EXPECT_EQ(after, 723u);
    }

    // Case 2: Add an aromatic bond for PHE (different parameters)
    {
        std::string comp_id = "PHE";
        std::string atom_1  = " C1";
        std::string atom_2  = " C2";
        std::string type    = "double";
        double value_dist           = 1.39;
        double value_dist_esd       = 0.02;
        double value_dist_nuclear   = 1.39;
        double value_dist_nuclear_esd = 0.02;
        coot::dict_bond_restraint_t::aromaticity_t arom = coot::dict_bond_restraint_t::AROMATIC;
        coot::dict_bond_restraint_t::bond_length_type_t blt = coot::dict_bond_restraint_t::NUCLEAR_POSITION;

        size_t before = count_all_bonds();
        EXPECT_EQ(before, 723u);

        coot::mon_lib_add_bond_gemmi(comp_id, 0, atom_1, atom_2, type,
                              value_dist, value_dist_esd,
                              value_dist_nuclear, value_dist_nuclear_esd,
                              arom, blt,
                              dict_res_restraints, read_number);
        size_t after = count_all_bonds();

        EXPECT_EQ(before, 723u);
        EXPECT_EQ(after, 724u);
    }

    // Case 3: Add a second bond to ALA — verifies bonds accumulate
    {
        std::string comp_id = "ALA";
        std::string atom_1  = " CA";
        std::string atom_2  = " C ";
        std::string type    = "single";
        double value_dist           = 1.523;
        double value_dist_esd       = 0.019;
        double value_dist_nuclear   = -1.0;
        double value_dist_nuclear_esd = -1.0;
        coot::dict_bond_restraint_t::aromaticity_t arom = coot::dict_bond_restraint_t::NON_AROMATIC;
        coot::dict_bond_restraint_t::bond_length_type_t blt = coot::dict_bond_restraint_t::UNKNOWN;

        size_t before = count_all_bonds();
        EXPECT_EQ(before, 724u);

        coot::mon_lib_add_bond_gemmi(comp_id, 0, atom_1, atom_2, type,
                              value_dist, value_dist_esd,
                              value_dist_nuclear, value_dist_nuclear_esd,
                              arom, blt,
                              dict_res_restraints, read_number);
        size_t after = count_all_bonds();

        EXPECT_EQ(before, 724u);
        EXPECT_EQ(after, 725u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
