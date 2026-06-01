#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, correct_link_distances) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Case 1: molecule with no links - function should be a no-op
    {
        EXPECT_EQ(st.models.size(), 1u);
        std::size_t n_links_before = st.connections.size();
        EXPECT_EQ(n_links_before, 0u);

        coot::util::correct_link_distances_gemmi(st);

        std::size_t n_links_after = st.connections.size();
        EXPECT_EQ(n_links_after, 0u);
    }

    // Case 2: add a link with deliberately wrong distance, then correct it
    {
        gemmi::Model& model = st.models[0];

        struct AtomInfo {
            gemmi::Atom* atom = nullptr;
            std::string chain_name;
            gemmi::ResidueId resid;
        };
        AtomInfo info1, info2;

        for (auto& chain : model.chains) {
            if (chain.name.empty()) continue;
            for (auto& res : chain.residues) {
                if (res.atoms.empty()) continue;
                if (!info1.atom) {
                    info1 = {&res.atoms[0], chain.name, res.seqid};
                } else if (!info2.atom) {
                    if (res.atoms.size() > 1) {
                        info2 = {&res.atoms[1], chain.name, res.seqid};
                    } else {
                        info2 = {&res.atoms[0], chain.name, res.seqid};
                    }
                }
            }
        }

        gemmi::Atom* a1 = info1.atom;
        gemmi::Atom* a2 = info2.atom;
        if (!a1 || !a2) {
            return;
        }

        gemmi::Vec3 diff = a1->pos - a2->pos;
        double real_dist = diff.length();

        // Build connection with wrong distance
        gemmi::Connection conn;
        conn.reported_distance = real_dist + 5.0;
        conn.type = gemmi::Connection::Type::Covale;

        conn.partner1 = gemmi::AtomAddress(info1.chain_name, info1.resid, a1->name, a1->altloc);
        conn.partner2 = gemmi::AtomAddress(info2.chain_name, info2.resid, a2->name, a2->altloc);

        st.connections.push_back(conn);
        EXPECT_NEAR(st.connections.back().reported_distance, real_dist + 5.0, 1e-4);

        coot::util::correct_link_distances_gemmi(st);

        EXPECT_NEAR(st.connections.back().reported_distance, real_dist, 0.01);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}