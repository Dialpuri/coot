#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// ... TEST(...) blocks here ...

TEST(OracleTest, add_guanidinium_hydrogens) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    gemmi::Model& model = st.models[0];

    ASSERT_FALSE(model.chains.empty());

    gemmi::Chain& chain = model.chains[0];

    // ARG 126
    {
        gemmi::Residue* res = nullptr;
        for (gemmi::Residue& r : chain.residues) {
            if (r.seqid.num.value == 126 && r.name == "ARG") {
                res = &r;
                break;
            }
        }
        ASSERT_NE(res, nullptr);

        std::cout << "ARG 126 atoms (" << res->atoms.size() << "):" << std::endl;
        for (int i = 0; i < static_cast<int>(res->atoms.size()); i++) {
            gemmi::Atom& a = res->atoms[i];
            std::cout << "  " << i << ": '" << a.name << "' x=" << a.pos.x << " y=" << a.pos.y << " z=" << a.pos.z << std::endl;
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
