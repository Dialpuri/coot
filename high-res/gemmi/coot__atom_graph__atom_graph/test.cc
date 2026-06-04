#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_graph) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Get atoms from the first residue of chain A to build minimal test data
    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];
    gemmi::Residue& res = chain.residues[0];

    EXPECT_STREQ(chain.name.c_str(), "A");
    EXPECT_STREQ(res.name.c_str(), "MET");
    EXPECT_EQ(res.seqid.num.value, 1);

    int n_atoms_in_res = res.atoms.size();
    EXPECT_EQ(n_atoms_in_res, 8);

    // Build coords from first few atoms (up to n_atoms_in_res, max 5)
    int n_use = (n_atoms_in_res < 5) ? n_atoms_in_res : 5;

    std::vector<gemmi::Vec3> coords;
    for (int i = 0; i < n_use; i++) {
        gemmi::Atom& at = res.atoms[i];
        coords.push_back(gemmi::Vec3(at.pos.x, at.pos.y, at.pos.z));
    }

    // Build connection_indices: each atom connected to the next (chain topology)
    std::vector<std::vector<coot::node_info>> connection_indices(n_use);
    for (int i = 0; i < n_use; i++) {
        if (i + 1 < n_use) {
            connection_indices[i].push_back(coot::node_info(i + 1));
        }
        if (i > 0) {
            connection_indices[i].push_back(coot::node_info(i - 1));
        }
    }

    EXPECT_EQ(connection_indices.size(), 5);
    EXPECT_EQ(coords.size(), 5);

    // Construct atom_graph
    coot::atom_graph ag;
    coot::atom_graph::atom_graph_gemmi(&st, connection_indices, coords, ag);

    EXPECT_EQ(ag.atom_info.size(), 5);
    EXPECT_EQ(ag.nodes.size(), 5);
    EXPECT_EQ(ag.spgr, std::string("P 21 21 21"));
    EXPECT_NEAR(ag.cell[0], 53.512, 1e-4);
    EXPECT_NEAR(ag.cell[1], 71.311, 1e-4);
    EXPECT_NEAR(ag.cell[2], 72.101, 1e-4);

    // Edge case: mismatched sizes (should produce error and empty atom_info/nodes)
    {
        std::vector<std::vector<coot::node_info>> bad_conn(1);
        std::vector<gemmi::Vec3> bad_coords(3);

        EXPECT_EQ(bad_conn.size(), 1);
        EXPECT_EQ(bad_coords.size(), 3);

        coot::atom_graph ag_bad;
        coot::atom_graph::atom_graph_gemmi(&st, bad_conn, bad_coords, ag_bad);
        EXPECT_EQ(ag_bad.atom_info.size(), 0);
        EXPECT_EQ(ag_bad.nodes.size(), 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}