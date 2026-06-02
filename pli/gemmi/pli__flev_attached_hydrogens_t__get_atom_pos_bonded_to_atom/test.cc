#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <cmath>
#include <string>
#include <vector>
#include "function.hh"

// Helper to trim atom name (remove leading/trailing spaces)
static std::string trim_atom_name(std::string name) {
    size_t start = name.find_first_not_of(' ');
    size_t end = name.find_last_not_of(' ');
    if (start == std::string::npos) return "";
    return name.substr(start, end - start + 1);
}

// Helper to find atom by trimmed name in a gemmi residue
static const gemmi::Atom* find_atom(const gemmi::Residue& res, const std::string& needle) {
    for (const gemmi::Atom& atom : res.atoms) {
        if (trim_atom_name(atom.name) == needle) {
            return &atom;
        }
    }
    return nullptr;
}

TEST(OracleTest, get_atom_pos_bonded_to_atom) {
    // Setup: load PDB and initialize geometry
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    coot::protein_geometry geom;
    geom.init_standard();

    // Find ILE residue //A/10
    const gemmi::Chain* chain_a = nullptr;
    for (const auto& model : st.models)
        for (const auto& chain : model.chains)
            if (chain.name == "A") {
                chain_a = &chain;
                break;
            }
    ASSERT_NE(chain_a, nullptr);

    const gemmi::Residue* ile_res = nullptr;
    for (const auto& res : chain_a->residues) {
        if (res.seqid.num.value == 10 && res.name == "ILE") {
            ile_res = &res;
            break;
        }
    }
    ASSERT_NE(ile_res, nullptr);
    EXPECT_EQ(ile_res->name, "ILE");

    // Create dictionary restraints for ILE
    std::pair<bool, coot::dictionary_residue_restraints_t> dict_pair =
        geom.get_monomer_restraints_at_least_minimal("ILE", 0);
    EXPECT_TRUE(dict_pair.first);

    // Create flev_attached_hydrogens_t instance
    pli::flev_attached_hydrogens_t flev(dict_pair.second);

    // Case 1 (oracle ground truth): CB bonded to HB in ILE
    // Function finds atom bonded to CB that is NOT HB (e.g. CA, CG1, CG2)
    {
        const gemmi::Atom* cb = find_atom(*ile_res, "CB");
        const gemmi::Atom* hb = find_atom(*ile_res, "HB");
        ASSERT_NE(cb, nullptr);
        ASSERT_NE(hb, nullptr);

        clipper::Coord_orth result = flev.get_atom_pos_bonded_to_atom_gemmi(cb, hb, ile_res, geom);

        EXPECT_NEAR(result.x(), 27.175, 1e-3);
        EXPECT_NEAR(result.y(), 12.237, 1e-3);
        EXPECT_NEAR(result.z(), 67.71, 1e-3);
    }

    // Case 2 (complementary): CA bonded to H in ILE — different bonded pair, different branch
    {
        const gemmi::Atom* ca = find_atom(*ile_res, "CA");
        const gemmi::Atom* h = find_atom(*ile_res, "H");
        if (ca && h) {
            clipper::Coord_orth result = flev.get_atom_pos_bonded_to_atom_gemmi(ca, h, ile_res, geom);
            EXPECT_FALSE(std::isnan(result.x()));
            EXPECT_FALSE(std::isnan(result.y()));
            EXPECT_FALSE(std::isnan(result.z()));
            EXPECT_TRUE(std::isfinite(result.x()));
            EXPECT_TRUE(std::isfinite(result.y()));
            EXPECT_TRUE(std::isfinite(result.z()));
        }
    }

    // Case 3 (complementary): CG1 bonded to HG1 in ILE — tests methyl-type hydrogen
    {
        const gemmi::Atom* cg1 = find_atom(*ile_res, "CG1");
        const gemmi::Atom* hg1 = find_atom(*ile_res, "HG1");
        if (cg1 && hg1) {
            clipper::Coord_orth result = flev.get_atom_pos_bonded_to_atom_gemmi(cg1, hg1, ile_res, geom);
            EXPECT_FALSE(std::isnan(result.x()));
            EXPECT_FALSE(std::isnan(result.y()));
            EXPECT_FALSE(std::isnan(result.z()));
            EXPECT_TRUE(std::isfinite(result.x()));
            EXPECT_TRUE(std::isfinite(result.y()));
            EXPECT_TRUE(std::isfinite(result.z()));
        }
    }

    // Case 4 (complementary): CG2 bonded to HG2 in ILE — tests another methyl hydrogen
    {
        const gemmi::Atom* cg2 = find_atom(*ile_res, "CG2");
        const gemmi::Atom* hg2 = find_atom(*ile_res, "HG2");
        if (cg2 && hg2) {
            clipper::Coord_orth result = flev.get_atom_pos_bonded_to_atom_gemmi(cg2, hg2, ile_res, geom);
            EXPECT_FALSE(std::isnan(result.x()));
            EXPECT_FALSE(std::isnan(result.y()));
            EXPECT_FALSE(std::isnan(result.z()));
            EXPECT_TRUE(std::isfinite(result.x()));
            EXPECT_TRUE(std::isfinite(result.y()));
            EXPECT_TRUE(std::isfinite(result.z()));
        }
    }

    // Case 5 (edge): HOH residue not in dictionary → should throw runtime_error
    {
        const gemmi::Chain* chain_b = nullptr;
        for (const auto& model : st.models)
            for (const auto& chain : model.chains)
                if (chain.name == "B") {
                    chain_b = &chain;
                    break;
                }
        if (chain_b) {
            const gemmi::Residue* hoh_res = nullptr;
            for (const auto& res : chain_b->residues) {
                if (res.seqid.num.value == 5001) {
                    hoh_res = &res;
                    break;
                }
            }
            if (hoh_res && !hoh_res->atoms.empty()) {
                EXPECT_THROW(
                    flev.get_atom_pos_bonded_to_atom_gemmi(&hoh_res->atoms[0], &hoh_res->atoms[0], hoh_res, geom),
                    std::runtime_error);
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}