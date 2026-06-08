#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mon_lib_add_angle) {
    // Standalone data structures mirroring protein_geometry internals
    std::vector<std::pair<int, coot::dictionary_residue_restraints_t>> dict_res_restraints;
    int read_number = 0;

    // Seed 41 entries with various comp_ids (none "ALA") to simulate init_standard()
    {
        const std::string nineteen[] = {
            "ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        };
        // 38 entries: 19 × 2 reps, plus 3 extra = 41
        for (int rep = 0; rep < 2; rep++) {
            for (unsigned int i = 0; i < 19; i++) {
                coot::dictionary_residue_restraints_t rest(nineteen[i], read_number);
                dict_res_restraints.push_back(std::make_pair(coot::protein_geometry::IMOL_ENC_AUTO, rest));
                read_number++;
            }
        }
        // Three more entries to reach 41
        coot::dictionary_residue_restraints_t rest1("WAT", read_number);
        dict_res_restraints.push_back(std::make_pair(coot::protein_geometry::IMOL_ENC_AUTO, rest1));
        read_number++;
        coot::dictionary_residue_restraints_t rest2("HOH", read_number);
        dict_res_restraints.push_back(std::make_pair(coot::protein_geometry::IMOL_ENC_AUTO, rest2));
        read_number++;
        coot::dictionary_residue_restraints_t rest3("MOL", read_number);
        dict_res_restraints.push_back(std::make_pair(coot::protein_geometry::IMOL_ENC_AUTO, rest3));
        read_number++;
    }

    // Case 1: add a standard angle restraint for ALA (N-CA-C angle)
    {
        std::string comp_id = "ALA";
        int imol_enc = coot::protein_geometry::IMOL_ENC_AUTO;
        std::string atom_id_1 = "N";
        std::string atom_id_2 = "CA";
        std::string atom_id_3 = "C";
        double value_angle = 110.0;
        double value_angle_esd = 2.0;

        size_t before = dict_res_restraints.size();
        EXPECT_EQ(before, 41u);

        coot::mon_lib_add_angle_gemmi(comp_id, imol_enc, atom_id_1, atom_id_2, atom_id_3,
                               value_angle, value_angle_esd,
                               dict_res_restraints, read_number);

        size_t after = dict_res_restraints.size();
        EXPECT_EQ(after, 42u);
        EXPECT_TRUE(after > before);
    }

    // Case 2: add a different angle for GLY (CA-C-O angle)
    {
        std::string comp_id = "GLY";
        int imol_enc = 1;
        std::string atom_id_1 = "CA";
        std::string atom_id_2 = "C";
        std::string atom_id_3 = "O";
        double value_angle = 120.5;
        double value_angle_esd = 1.5;

        size_t before = dict_res_restraints.size();
        EXPECT_EQ(before, 42u);

        coot::mon_lib_add_angle_gemmi(comp_id, imol_enc, atom_id_1, atom_id_2, atom_id_3,
                               value_angle, value_angle_esd,
                               dict_res_restraints, read_number);

        size_t after = dict_res_restraints.size();
        EXPECT_EQ(after, 43u);
        EXPECT_TRUE(after > before);
    }

    // Case 3: edge case - unknown residue, zero angle
    {
        std::string comp_id = "UNX";
        int imol_enc = coot::protein_geometry::IMOL_ENC_AUTO;
        std::string atom_id_1 = "X1";
        std::string atom_id_2 = "X2";
        std::string atom_id_3 = "X3";
        double value_angle = 0.0;
        double value_angle_esd = 0.0;

        size_t before = dict_res_restraints.size();
        EXPECT_EQ(before, 43u);

        coot::mon_lib_add_angle_gemmi(comp_id, imol_enc, atom_id_1, atom_id_2, atom_id_3,
                               value_angle, value_angle_esd,
                               dict_res_restraints, read_number);

        size_t after = dict_res_restraints.size();
        EXPECT_EQ(after, 44u);
        EXPECT_TRUE(after > before);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
