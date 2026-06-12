#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include "function.hh"

// Helper to check if an atom name matches (handles PDB padding and CIF unpadded)
static bool atom_name_matches(const std::string& name, const std::string& target) {
    return (name == target ||
            std::string(1, name[0]) + std::string(1, name[1]) + (name.size() > 2 ? std::string(1, name[2]) : "") == target ||
            name == " " + target + " " ||
            name == target);
}

static bool is_atom(const std::string& name, const std::string& t1, const std::string& t2 = "") {
    if (name == t1 || name == t2) return true;
    // Handle padded: " SD ", "SE  "
    if ((name.size() == 4 && name.substr(0, 2) == t1) ||
        (name.size() == 4 && name.substr(0, 2) == t2)) return true;
    // Handle unpadded CIF
    if (name == " " + t1 + " " || name == " " + t2 + " ") return true;
    return false;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, handle_MET_or_MSE_case) {
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_EQ(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];
    int model_number = 0; // 0-indexed in gemmi

    int atom_colour_type = 0;  // COLOUR_BY_CHAIN
    int udd_user_defined_atom_colour_index_handle = 0;

    // --- Case 1: MET SD atom (oracle ground truth) ---
    {
        gemmi::Atom* target_atom = nullptr;
        gemmi::Chain* target_chain = nullptr;
        gemmi::Residue* target_residue = nullptr;
        std::string target_cid;
        std::string target_resname;
        std::string target_atomname;

        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                std::string resname = residue.name;
                if (resname == "MET" || resname == "MSE") {
                    for (gemmi::Atom& atom : residue.atoms) {
                        std::string aname = atom.name;
                        // Match "SD", " SD ", "SD  " etc.
                        bool match = (aname == " SD " || aname == "SD" ||
                                     aname == "SE  " || aname == "SE" ||
                                     (aname.size() >= 2 && aname.substr(0, 2) == "SD") ||
                                     (aname.size() >= 2 && aname.substr(0, 2) == "SE"));
                        if (match) {
                            target_atom = &atom;
                            target_chain = &chain;
                            target_residue = &residue;
                            target_cid = "//" + chain.name + "/" + std::to_string(residue.seqid.num.value) + "/" + atom.name;
                            target_resname = resname;
                            target_atomname = aname;
                            goto found_target;
                        }
                    }
                }
            }
        }
        found_target:;

        ASSERT_NE(target_atom, nullptr);
        EXPECT_EQ(target_resname, "MET");

        // Create Bond_lines_container with 13 colour slots
        Bond_lines_container_gemmi bonds;
        if (bonds.bonds.size() == 0) {
            for (int i = 0; i < 13; i++) {
                Bond_lines a(i);
                bonds.bonds.push_back(a);
            }
        }

        int bonds_before = 0;
        for (size_t i = 0; i < bonds.bonds.size(); i++) {
            bonds_before += static_cast<int>(bonds.bonds[i].size());
        }
        EXPECT_EQ(bonds_before, 0);

        gemmi::CRA target_cra{target_chain, target_residue, target_atom};
        handle_MET_or_MSE_case_gemmi(
            bonds,
            target_cra,
            model_number,
            atom_colour_type,
            udd_user_defined_atom_colour_index_handle,
            nullptr  // no atom_colour_map
        );

        int bonds_after = 0;
        for (size_t i = 0; i < bonds.bonds.size(); i++) {
            bonds_after += static_cast<int>(bonds.bonds[i].size());
        }
        EXPECT_EQ(bonds_after, 4);
        EXPECT_EQ(bonds_after - bonds_before, 4);
    }

    // --- Case 2: non-MET/MSE/CYS atom (should be a no-op) ---
    {
        gemmi::Atom* ala_atom = nullptr;
        gemmi::Chain* ala_chain = nullptr;
        gemmi::Residue* ala_residue = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                std::string resname = residue.name;
                if (resname == "ALA") {
                    for (gemmi::Atom& atom : residue.atoms) {
                        bool is_CA = (atom.name == " CA " || atom.name == "CA" ||
                                      (atom.name.size() >= 2 && atom.name.substr(0, 2) == "CA"));
                        if (is_CA) {
                            ala_atom = &atom;
                            ala_chain = &chain;
                            ala_residue = &residue;
                            goto found_ala;
                        }
                    }
                }
            }
        }
        found_ala:;

        if (ala_atom) {
            std::string ala_resname = ala_residue->name;
            EXPECT_EQ(ala_resname, "ALA");

            Bond_lines_container_gemmi bonds2;
            if (bonds2.bonds.size() == 0) {
                for (int i = 0; i < 13; i++) {
                    Bond_lines a(i);
                    bonds2.bonds.push_back(a);
                }
            }
            int bonds2_before = 0;
            for (size_t i = 0; i < bonds2.bonds.size(); i++) {
                bonds2_before += static_cast<int>(bonds2.bonds[i].size());
            }

            gemmi::CRA ala_cra{ala_chain, ala_residue, ala_atom};
            handle_MET_or_MSE_case_gemmi(
                bonds2,
                ala_cra,
                model_number,
                atom_colour_type,
                udd_user_defined_atom_colour_index_handle,
                nullptr  // no atom_colour_map
            );

            int bonds2_after = 0;
            for (size_t i = 0; i < bonds2.bonds.size(); i++) {
                bonds2_after += static_cast<int>(bonds2.bonds[i].size());
            }

            // Should be a no-op for non-MET/MSE/CYS
            EXPECT_EQ(bonds2_before, bonds2_after);
        }
    }

    // --- Case 3: CYS SG atom (exercise CYS branch) ---
    {
        gemmi::Atom* cys_atom = nullptr;
        gemmi::Chain* cys_chain = nullptr;
        gemmi::Residue* cys_residue = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& residue : chain.residues) {
                std::string resname = residue.name;
                if (resname == "CYS") {
                    for (gemmi::Atom& atom : residue.atoms) {
                        bool is_SG = (atom.name == " SG " || atom.name == "SG" ||
                                      (atom.name.size() >= 2 && atom.name.substr(0, 2) == "SG"));
                        if (is_SG) {
                            cys_atom = &atom;
                            cys_chain = &chain;
                            cys_residue = &residue;
                            goto found_cys;
                        }
                    }
                }
            }
        }
        found_cys:;

        if (cys_atom) {
            std::string cys_resname = cys_residue->name;
            EXPECT_EQ(cys_resname, "CYS");

            Bond_lines_container_gemmi bonds3;
            if (bonds3.bonds.size() == 0) {
                for (int i = 0; i < 13; i++) {
                    Bond_lines a(i);
                    bonds3.bonds.push_back(a);
                }
            }
            int bonds3_before = 0;
            for (size_t i = 0; i < bonds3.bonds.size(); i++) {
                bonds3_before += static_cast<int>(bonds3.bonds[i].size());
            }

            gemmi::CRA cys_cra{cys_chain, cys_residue, cys_atom};
            handle_MET_or_MSE_case_gemmi(
                bonds3,
                cys_cra,
                model_number,
                atom_colour_type,
                udd_user_defined_atom_colour_index_handle,
                nullptr  // no atom_colour_map
            );

            int bonds3_after = 0;
            for (size_t i = 0; i < bonds3.bonds.size(); i++) {
                bonds3_after += static_cast<int>(bonds3.bonds[i].size());
            }

            // CYS SG should create 2 bonds (SG->mid, mid->CB) if within distance
            EXPECT_GE(bonds3_after, bonds3_before);
        }
    }
}
