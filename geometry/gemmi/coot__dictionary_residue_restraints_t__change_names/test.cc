#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <cctype>
#include <string>
#include <vector>
#include <utility>
#include "function.hh"

static std::string pad_atom_name_pdb_test(const std::string& name) {
    if (name.size() >= 4) return name.substr(0, 4);
    if (!name.empty() && std::isalpha(static_cast<unsigned char>(name[0]))) {
        std::string padded = " " + name;
        while (padded.size() < 4) padded += " ";
        return padded;
    }
    std::string padded = name;
    while (padded.size() < 4) padded += " ";
    return padded;
}

TEST(OracleTest, dictionary_residue_restraints_t_change_names) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: valid rename — change N -> NX on MET residue at //A/1
    {
        gemmi::Residue* res = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& r : chain.residues) {
                    if (r.seqid.num.value == 1) {
                        res = &r;
                        break;
                    }
                }
                break;
            }
        }
        ASSERT_NE(res, nullptr);

        coot::dictionary_residue_restraints_t dict_restraints;

        // Capture atom names BEFORE
        std::string atoms_before;
        for (size_t i = 0; i < res->atoms.size(); i++) {
            if (i > 0) atoms_before += ",";
            atoms_before += pad_atom_name_pdb_test(res->atoms[i].name);
        }
        std::string resname_before = res->name;

        EXPECT_EQ(resname_before, "MET");
        EXPECT_EQ(atoms_before, " N  , CA , C  , O  , CB , CG , SD , CE ");

        std::vector<std::pair<std::string, std::string>> changes;
        changes.push_back(std::make_pair("N", "NX"));
        std::string new_comp_id = "XAA";

        bool changed = dict_restraints.change_names_gemmi(res, changes, new_comp_id);

        // Capture after
        std::string atoms_after;
        for (size_t i = 0; i < res->atoms.size(); i++) {
            if (i > 0) atoms_after += ",";
            atoms_after += pad_atom_name_pdb_test(res->atoms[i].name);
        }
        std::string resname_after = res->name;

        EXPECT_FALSE(changed);
        EXPECT_EQ(resname_after, "MET");
        EXPECT_EQ(atoms_after, " N  , CA , C  , O  , CB , CG , SD , CE ");
    }

    // Case 2: empty changes — should not change anything
    {
        gemmi::Residue* res2 = nullptr;
        for (auto& chain : model.chains) {
            if (chain.name == "A") {
                for (auto& r : chain.residues) {
                    if (r.seqid.num.value == 2) {
                        res2 = &r;
                        break;
                    }
                }
                break;
            }
        }
        ASSERT_NE(res2, nullptr);

        std::string atoms2_before;
        for (size_t i = 0; i < res2->atoms.size(); i++) {
            if (i > 0) atoms2_before += ",";
            atoms2_before += pad_atom_name_pdb_test(res2->atoms[i].name);
        }
        std::string resname2_before = res2->name;

        EXPECT_EQ(resname2_before, "GLU");
        EXPECT_EQ(atoms2_before, " N  , CA , C  , O  , CB , CG , CD , OE1, OE2");

        coot::dictionary_residue_restraints_t dict_restraints;
        std::vector<std::pair<std::string, std::string>> empty_changes;
        std::string new_comp_id2 = "YYY";

        bool changed2 = dict_restraints.change_names_gemmi(res2, empty_changes, new_comp_id2);

        std::string atoms2_after;
        for (size_t i = 0; i < res2->atoms.size(); i++) {
            if (i > 0) atoms2_after += ",";
            atoms2_after += pad_atom_name_pdb_test(res2->atoms[i].name);
        }
        std::string resname2_after = res2->name;

        EXPECT_FALSE(changed2);
        EXPECT_EQ(resname2_after, "GLU");
        EXPECT_EQ(atoms2_after, " N  , CA , C  , O  , CB , CG , CD , OE1, OE2");
    }

    // Case 3: null residue — should return false without crashing
    {
        coot::dictionary_residue_restraints_t dict_restraints;
        std::vector<std::pair<std::string, std::string>> changes3;
        changes3.push_back(std::make_pair("CA", "CB"));
        std::string new_comp_id3 = "ZZZ";

        bool changed3 = dict_restraints.change_names_gemmi(nullptr, changes3, new_comp_id3);

        EXPECT_FALSE(changed3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}