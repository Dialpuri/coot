#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, best_fit_phi_psi_attaching_oxygen_position_update) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model &model = st.models[0];

    // case: oracle case 1 — standard residue with manually added atoms at (0,0,0)
    {
        const gemmi::Residue *residue_p = nullptr;
        for (gemmi::Chain &chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        residue_p = &res;
                        break;
                    }
                }
            }
        }
        ASSERT_NE(residue_p, nullptr);

        // Find chain and residue info
        const gemmi::Chain *chain_p = nullptr;
        for (gemmi::Chain &chain : model.chains) {
            if (chain.name == "A") {
                chain_p = &chain;
                break;
            }
        }
        ASSERT_NE(chain_p, nullptr);

        std::string chain_id = chain_p->name;
        int resno = residue_p->seqid.num.value;
        EXPECT_EQ(chain_id, "A");
        EXPECT_EQ(resno, 10);
        EXPECT_STREQ(residue_p->name.c_str(), "ILE");

        // Get CA and C positions
        clipper::Coord_orth ca_pos(0, 0, 0), c_pos(0, 0, 0);
        bool found_ca = false, found_c = false;
        for (const gemmi::Atom &atom : residue_p->atoms) {
            if (atom.name == " CA ") {
                ca_pos = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
                found_ca = true;
            }
            if (atom.name == " C  ") {
                c_pos = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
                found_c = true;
            }
        }
        EXPECT_FALSE(found_ca);
        EXPECT_FALSE(found_c);
        EXPECT_NEAR(ca_pos.x(), 0, 1e-4);
        EXPECT_NEAR(ca_pos.y(), 0, 1e-4);
        EXPECT_NEAR(ca_pos.z(), 0, 1e-4);
        EXPECT_NEAR(c_pos.x(), 0, 1e-4);
        EXPECT_NEAR(c_pos.y(), 0, 1e-4);
        EXPECT_NEAR(c_pos.z(), 0, 1e-4);

        // Build mmol: fragment with 1 residue, atoms at [0]=CA, [1]=C, [2]=N
        coot::minimol::fragment frag;
        frag.residues.push_back(coot::minimol::residue(0, "ALA"));
        frag.residues[0].addatom(" CA ", "C", ca_pos, "", 20, 1);
        frag.residues[0].addatom(" C  ", "C", c_pos, "", 20, 1);
        frag.residues[0].addatom(" N  ", "N", clipper::Coord_orth(10, 20, 30), "", 20, 1);

        coot::minimol::molecule mm;
        mm.fragments.push_back(frag);

        EXPECT_EQ(mm.fragments.size(), 1u);
        EXPECT_EQ(mm.fragments[0].residues[0].atoms.size(), 3u);

        clipper::Coord_orth pos = coot::residue_by_phi_psi::best_fit_phi_psi_attaching_oxygen_position_update_gemmi(mm, residue_p);

        EXPECT_NEAR(pos.x(), 0, 1e-4);
        EXPECT_NEAR(pos.y(), 0, 1e-4);
        EXPECT_NEAR(pos.z(), 0, 1e-4);
    }

    // case: oracle case 2 — zero-fragment molecule
    {
        const gemmi::Residue *residue_p = nullptr;
        for (gemmi::Chain &chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        residue_p = &res;
                        break;
                    }
                }
            }
        }
        ASSERT_NE(residue_p, nullptr);

        coot::minimol::molecule mm_empty;
        EXPECT_EQ(mm_empty.fragments.size(), 0u);

        clipper::Coord_orth pos = coot::residue_by_phi_psi::best_fit_phi_psi_attaching_oxygen_position_update_gemmi(mm_empty, residue_p);

        EXPECT_NEAR(pos.x(), 0, 1e-4);
        EXPECT_NEAR(pos.y(), 0, 1e-4);
        EXPECT_NEAR(pos.z(), 0, 1e-4);
    }

    // case: complementary — real coords extracted from PDB residue, should exercise core logic
    {
        const gemmi::Residue *residue_p = nullptr;
        for (gemmi::Chain &chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue &res : chain.residues) {
                    if (res.seqid.num.value == 10) {
                        residue_p = &res;
                        break;
                    }
                }
            }
        }
        ASSERT_NE(residue_p, nullptr);

        // Extract CA, C, N positions from the real residue
        EXPECT_GT(residue_p->atoms.size(), 0u);

        clipper::Coord_orth ca_pos(0, 0, 0), c_pos(0, 0, 0), n_pos(0, 0, 0);
        bool has_ca = false, has_c = false, has_n = false;
        for (const gemmi::Atom &atom : residue_p->atoms) {
            std::string aname = atom.name;
            if (aname == "CA" || aname == " CA ") {
                ca_pos = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
                has_ca = true;
            }
            if (aname == "C" || aname == " C  ") {
                c_pos = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
                has_c = true;
            }
            if (aname == "N" || aname == " N  ") {
                n_pos = clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z);
                has_n = true;
            }
        }
        EXPECT_TRUE(has_ca);
        EXPECT_TRUE(has_c);
        EXPECT_TRUE(has_n);

        // Build minimol with real coords: atoms[0]=CA, atoms[1]=C, atoms[2]=N
        coot::minimol::fragment frag;
        frag.residues.push_back(coot::minimol::residue(0, "ILE"));
        frag.residues[0].addatom("CA", "C", ca_pos, "", 1, 1);
        frag.residues[0].addatom("C", "C", c_pos, "", 1, 1);
        frag.residues[0].addatom("N", "N", n_pos, "", 1, 1);

        coot::minimol::molecule mm;
        mm.fragments.push_back(frag);

        EXPECT_EQ(mm.fragments.size(), 1u);
        EXPECT_EQ(mm.fragments[0].residues[0].atoms.size(), 3u);

        clipper::Coord_orth pos = coot::residue_by_phi_psi::best_fit_phi_psi_attaching_oxygen_position_update_gemmi(mm, residue_p);

        // With real CA/C/N coords the function should compute a non-trivial O position
        EXPECT_FALSE(std::isnan(pos.x()));
        EXPECT_FALSE(std::isnan(pos.y()));
        EXPECT_FALSE(std::isnan(pos.z()));
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}