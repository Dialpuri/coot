#include <gtest/gtest.h>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

// Helper: find a residue by chain ID string and sequence number in a gemmi Structure
static gemmi::Residue const *find_residue(gemmi::Structure const &st, std::string chain_id, int seq_num) {
    for (gemmi::Model const &model : st.models) {
        for (gemmi::Chain const &chain : model.chains) {
            if (chain.name == chain_id) {
                for (gemmi::Residue const &res : chain.residues) {
                    if (res.seqid.num.value == seq_num)
                        return &res;
                }
            }
        }
    }
    return nullptr;
}

// Helper: pad atom name to 4 chars (left-padded with spaces, MMDB style)
static std::string atom_name_4c(std::string const &name) {
    std::string trimmed = name;
    size_t a = trimmed.find_first_not_of(' ');
    size_t b = trimmed.find_last_not_of(' ');
    if (a == std::string::npos) return "    ";
    trimmed = trimmed.substr(a, b - a + 1);
    int ilen = static_cast<int>(trimmed.length());
    if (ilen >= 4) return trimmed.substr(0, 4);
    return std::string(4 - ilen, ' ') + trimmed;
}

// Helper: check if an atom is nitrogen (handles CIF element parsing edge cases)
static bool is_nitrogen(gemmi::Atom const &atom) {
    std::string ele = atom.element.name();
    // Trim trailing whitespace
    while (!ele.empty() && ele.back() == ' ') ele.pop_back();
    if (ele == "N") return true;
    // Fallback: atom name starts with 'N' (common for N atoms in non-standard residues)
    std::string trimmed = atom.name;
    size_t a = trimmed.find_first_not_of(' ');
    if (a != std::string::npos && trimmed[a] == 'N') return true;
    return false;
}

TEST(OracleTest, get_ligand_cations) {
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");

    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Ligand residue LZA at chain A/1299 — no library restraints → 0 cations
    {
        gemmi::Residue const *res = find_residue(st, "A", 1299);
        EXPECT_NE(res, nullptr);

        std::string res_name = res->name;
        EXPECT_EQ(res_name, "LZA");

        auto restraints_pair = geom.get_monomer_restraints(res_name, 0);
        coot::dictionary_residue_restraints_t monomer_restraints;
        if (restraints_pair.first) {
            monomer_restraints = restraints_pair.second;
        }
        EXPECT_EQ(monomer_restraints.bond_restraint.size(), 0u);

        std::vector<std::pair<std::string, clipper::Coord_orth>> cations =
            pli::pi_stacking_container_t::get_ligand_cations_gemmi(*res, monomer_restraints);

        EXPECT_EQ(cations.size(), 0u);
    }

    // Case 2: ARG residue — N atoms exist but none have 4 bonds → 0 cations
    {
        gemmi::Residue const *res = find_residue(st, "A", 50);
        if (!res) {
            res = find_residue(st, "A", 10);
        }
        EXPECT_NE(res, nullptr);

        std::string res_name = res->name;
        EXPECT_EQ(res_name, "ARG");

        auto restraints_pair = geom.get_monomer_restraints(res_name, 0);
        EXPECT_TRUE(restraints_pair.first);
        coot::dictionary_residue_restraints_t monomer_restraints = restraints_pair.second;
        EXPECT_EQ(monomer_restraints.bond_restraint.size(), 26u);

        std::vector<std::pair<std::string, clipper::Coord_orth>> cations =
            pli::pi_stacking_container_t::get_ligand_cations_gemmi(*res, monomer_restraints);

        EXPECT_EQ(cations.size(), 0u);
    }

    // Case 3: Synthesize restraints with N atom having 4 single bonds to real atoms → 1 cation
    {
        gemmi::Residue const *res = find_residue(st, "A", 1299);
        EXPECT_NE(res, nullptr);

        std::string res_name = res->name;
        EXPECT_EQ(res_name, "LZA");

        // Find N atoms and other atoms in LZA
        std::string nitrogen_name;
        std::vector<std::string> other_atom_names;
        for (gemmi::Atom const &atom : res->atoms) {
            if (is_nitrogen(atom)) {
                nitrogen_name = atom_name_4c(atom.name);
            } else {
                other_atom_names.push_back(atom_name_4c(atom.name));
            }
        }

        EXPECT_EQ(nitrogen_name, " N24");

        // Create synthetic restraints: N atom with 4 single bonds to existing atoms
        coot::dictionary_residue_restraints_t syn_restraints;
        for (int i = 0; i < 4; i++) {
            coot::dict_bond_restraint_t br(nitrogen_name,
                                           (i < static_cast<int>(other_atom_names.size()))
                                               ? other_atom_names[i]
                                               : other_atom_names[0],
                                           "single");
            syn_restraints.bond_restraint.push_back(br);
        }

        EXPECT_EQ(syn_restraints.bond_restraint.size(), 4u);

        std::vector<std::pair<std::string, clipper::Coord_orth>> cations =
            pli::pi_stacking_container_t::get_ligand_cations_gemmi(*res, syn_restraints);

        EXPECT_EQ(cations.size(), 1u);
        EXPECT_EQ(cations[0].first, " N24");
        EXPECT_NEAR(cations[0].second.x(), 32.96, 1e-1);
        EXPECT_NEAR(cations[0].second.y(), 9.451, 1e-1);
        EXPECT_NEAR(cations[0].second.z(), 66.12, 1e-1);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}