#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <string>
#include <vector>
#include <utility>
#include "clipper/core/coords.h"
#include "geometry/protein-geometry.hh"
#include "function.hh"

// Helper to get a residue by chain name and sequence number
static const gemmi::Residue* find_residue(
    const gemmi::Model& model, const std::string& chain_id, int seq_num) {
    for (const gemmi::Chain& chain : model.chains) {
        if (chain.name == chain_id) {
            for (const gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return &res;
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, update_chem_comp_atoms_from_residue) {
    {
        // Case 1: normal residue (MET) with populated atom_info
        gemmi::Structure st = gemmi::read_pdb_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        const gemmi::Model& model = st.models[0];

        const gemmi::Residue* residue_p = find_residue(model, "A", 1);
        EXPECT_NE(residue_p, nullptr);
        EXPECT_EQ(residue_p->name, "MET");

        // Count atoms in the residue
        size_t n_residue_atoms = residue_p->atoms.size();
        EXPECT_EQ(n_residue_atoms, 8u);

        // Create restraints object with atom_info entries for each non-terminal atom
        // In gemmi, there are no TER atoms, so all atoms are real atoms
        coot::dictionary_residue_restraints_t restraints;
        restraints.residue_info.comp_id = residue_p->name;

        for (size_t i = 0; i < n_residue_atoms; i++) {
            const gemmi::Atom& at = residue_p->atoms[i];
            std::string atom_name(at.name);

            coot::dict_atom da;
            da.atom_id_4c = atom_name;
            da.model_Cartn = std::make_pair(true, clipper::Coord_orth(-9999.0, -9999.0, -9999.0));
            restraints.atom_info.push_back(da);
        }

        EXPECT_EQ(restraints.atom_info.size(), 8u);

        // Record BEFORE coordinates for CA and C (carbonyl carbon)
        // Replicate MMDB matching: atom_id_4c starts with 'C' but not 'CA'/'CB'
        std::string ca_before_str = "NOT FOUND";
        std::string c_before_str = "NOT FOUND";
        for (size_t j = 0; j < restraints.atom_info.size(); j++) {
            const std::string& an = restraints.atom_info[j].atom_id_4c;
            if (an == "CA") {
                ca_before_str = restraints.atom_info[j].model_Cartn.second.format();
            }
            // Matches C, CE, CG, etc. (anything starting with 'C' but not CA/CB)
            if (!an.empty() && an[0] == 'C' && an != "CA" && an != "CB") {
                c_before_str = restraints.atom_info[j].model_Cartn.second.format();
            }
        }
        EXPECT_EQ(ca_before_str, "xyz = (     -9999,     -9999,     -9999)");
        EXPECT_EQ(c_before_str, "xyz = (     -9999,     -9999,     -9999)");

        // === CALL THE FUNCTION ===
        coot::update_chem_comp_atoms_from_residue_gemmi(*residue_p, &restraints);

        // Record AFTER coordinates
        std::string ca_after_str = "NOT FOUND";
        std::string c_after_str = "NOT FOUND";
        for (size_t j = 0; j < restraints.atom_info.size(); j++) {
            const std::string& an = restraints.atom_info[j].atom_id_4c;
            if (an == "CA") {
                ca_after_str = restraints.atom_info[j].model_Cartn.second.format();
            }
            // Matches C, CE, CG, etc. (anything starting with 'C' but not CA/CB)
            if (!an.empty() && an[0] == 'C' && an != "CA" && an != "CB") {
                c_after_str = restraints.atom_info[j].model_Cartn.second.format();
            }
        }

        // Oracle output values
        EXPECT_EQ(ca_after_str, "xyz = (      13.4,      3.78,      74.7)");
        EXPECT_EQ(c_after_str, "xyz = (     15.47,     7.548,     75.57)");
        EXPECT_NE(ca_before_str, ca_after_str);  // CA_updated: true
        EXPECT_NE(c_before_str, c_after_str);    // C_updated: true
    }

    {
        // Case 2: empty atom_info — no-op, nothing to update
        gemmi::Structure st = gemmi::read_pdb_file(
            "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        const gemmi::Model& model = st.models[0];

        const gemmi::Residue* residue_p = find_residue(model, "A", 1);
        EXPECT_NE(residue_p, nullptr);

        coot::dictionary_residue_restraints_t empty_restraints;
        EXPECT_EQ(empty_restraints.atom_info.size(), 0u);

        EXPECT_NO_THROW(
            coot::update_chem_comp_atoms_from_residue_gemmi(*residue_p, &empty_restraints));

        EXPECT_EQ(empty_restraints.atom_info.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}