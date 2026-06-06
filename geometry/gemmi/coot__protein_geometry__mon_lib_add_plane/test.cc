#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, mon_lib_add_plane) {
    // Standalone data structures mirroring protein_geometry internals
    std::vector<std::pair<int, coot::dictionary_residue_restraints_t>> dict_res_restraints;
    int read_number = 0;

    // Seed 41 entries with various comp_ids (none "ALA") to simulate init_standard()
    {
        const std::string nineteen[] = {
            "ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        };
        // 38 entries: 19 non-ALA × 2 reps, plus 3 extra non-ALA entries = 41
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

    // Helper to count plane restraints for a given comp_id
    auto count_planes = [&](const std::vector<std::pair<int, coot::dictionary_residue_restraints_t>> &drs, const std::string &comp_id, int imol_enc) -> unsigned int {
        unsigned int count = 0;
        for (unsigned int i = 0; i < drs.size(); i++) {
            if (drs[i].second.residue_info.comp_id == comp_id) {
                if (drs[i].first == imol_enc) {
                    count += drs[i].second.plane_restraint.size();
                }
            }
        }
        return count;
    };

    auto count_entries = [&]() -> unsigned int {
        return dict_res_restraints.size();
    };

    // Case 1: Add a plane to a new comp_id (not yet in dict_res_restraints)
    {
        std::string comp_id = "ALA";
        int imol_enc = coot::protein_geometry::IMOL_ENC_AUTO;
        std::string plane_id = "PL1";
        std::string atom_id = "CA";
        double dist_esd = 0.02;

        unsigned int before_entries = count_entries();
        unsigned int before_planes = count_planes(dict_res_restraints, comp_id, imol_enc);

        coot::mon_lib_add_plane_gemmi(comp_id, imol_enc, plane_id, atom_id, dist_esd, dict_res_restraints, read_number);

        unsigned int after_entries = count_entries();
        unsigned int after_planes = count_planes(dict_res_restraints, comp_id, imol_enc);

        EXPECT_EQ(before_entries, 41u);
        EXPECT_EQ(after_entries, 42u);
        EXPECT_EQ(before_planes, 0u);
        EXPECT_EQ(after_planes, 1u);
    }

    // Case 2: Add another atom to an existing plane (same comp_id, same plane_id)
    {
        std::string comp_id = "ALA";
        int imol_enc = coot::protein_geometry::IMOL_ENC_AUTO;
        std::string plane_id = "PL1";
        std::string atom_id = "C";
        double dist_esd = 0.015;

        unsigned int before_planes = count_planes(dict_res_restraints, comp_id, imol_enc);

        coot::mon_lib_add_plane_gemmi(comp_id, imol_enc, plane_id, atom_id, dist_esd, dict_res_restraints, read_number);

        unsigned int after_planes = count_planes(dict_res_restraints, comp_id, imol_enc);

        EXPECT_EQ(before_planes, 1u);
        EXPECT_EQ(after_planes, 1u);
    }

    // Case 3: Add a new plane_id to the same comp_id
    {
        std::string comp_id = "ALA";
        int imol_enc = coot::protein_geometry::IMOL_ENC_AUTO;
        std::string plane_id = "PL2";
        std::string atom_id = "N";
        double dist_esd = 0.01;

        unsigned int before_planes = count_planes(dict_res_restraints, comp_id, imol_enc);

        coot::mon_lib_add_plane_gemmi(comp_id, imol_enc, plane_id, atom_id, dist_esd, dict_res_restraints, read_number);

        unsigned int after_planes = count_planes(dict_res_restraints, comp_id, imol_enc);

        EXPECT_EQ(before_planes, 1u);
        EXPECT_EQ(after_planes, 2u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}