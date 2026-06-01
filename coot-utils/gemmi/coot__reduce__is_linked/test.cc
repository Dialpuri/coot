#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/metadata.hpp>
#include <gemmi/seqid.hpp>
#include "function.hh"

TEST(OracleTest, reduce_is_linked) {
    // Load a PDB with SSBOND links
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-dev/coot/data/tutorial-modern.pdb");

    gemmi::Model &model = st.models[0];

    // Find chain B
    gemmi::Chain *chain = nullptr;
    for (gemmi::Chain &ch : model.chains) {
        if (ch.name == "B") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    // Find CYS residues 7 and 96 (from SSBOND record)
    gemmi::Residue *cys7 = nullptr;
    gemmi::Residue *cys96 = nullptr;
    for (gemmi::Residue &res : chain->residues) {
        if (res.seqid.num.value == 7 && res.name == "CYS") {
            cys7 = &res;
        }
        if (res.seqid.num.value == 96 && res.name == "CYS") {
            cys96 = &res;
        }
    }

    ASSERT_NE(cys7, nullptr);
    ASSERT_NE(cys96, nullptr);

    // Find SG atoms (gemmi trims atom names, so " SG " becomes "SG")
    gemmi::Atom *sg7 = nullptr;
    gemmi::Atom *sg96 = nullptr;
    for (gemmi::Atom &atom : cys7->atoms) {
        if (atom.name == "SG") {
            sg7 = &atom;
            break;
        }
    }
    for (gemmi::Atom &atom : cys96->atoms) {
        if (atom.name == "SG") {
            sg96 = &atom;
            break;
        }
    }
    ASSERT_NE(sg7, nullptr);
    ASSERT_NE(sg96, nullptr);

    // Manually add a connection between SG atoms (SSBOND)
    gemmi::Connection conn;
    conn.type = gemmi::Connection::Disulf;

    // Build partner1 AtomAddress
    gemmi::ResidueId res_id7;
    res_id7.seqid = cys7->seqid;
    res_id7.name  = cys7->name;
    conn.partner1 = gemmi::AtomAddress(
        chain->name, res_id7, sg7->name, sg7->altloc);

    // Build partner2 AtomAddress
    gemmi::ResidueId res_id96;
    res_id96.seqid = cys96->seqid;
    res_id96.name  = cys96->name;
    conn.partner2 = gemmi::AtomAddress(
        chain->name, res_id96, sg96->name, sg96->altloc);

    st.connections.push_back(conn);

    // Case 1: SG atom of CYS B/7 — should be linked (true)
    {
        gemmi::CRA cra{chain, cys7, nullptr};
        std::string atom_name = " SG ";
        bool result = coot::reduce::is_linked_gemmi(atom_name, cra, st);
        EXPECT_TRUE(result);
    }

    // Case 2: CA atom of CYS B/7 — should NOT be linked (false)
    {
        gemmi::CRA cra{chain, cys7, nullptr};
        std::string atom_name = " CA ";
        bool result = coot::reduce::is_linked_gemmi(atom_name, cra, st);
        EXPECT_FALSE(result);
    }

    // Case 3: SG atom of CYS B/96 — should be linked (true, other end)
    {
        gemmi::CRA cra{chain, cys96, nullptr};
        std::string atom_name = " SG ";
        bool result = coot::reduce::is_linked_gemmi(atom_name, cra, st);
        EXPECT_TRUE(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}