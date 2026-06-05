#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "api/molecules-container.hh"
#include "clipper/core/xmap.h"
#include "clipper/core/coords.h"
#include <tuple>
#include <vector>
#include <string>

// Helper: find residue by chain ID and sequence number in a gemmi structure
static const gemmi::Residue* find_residue(const gemmi::Model& model,
                                           const std::string& chain_id,
                                           int seqnum) {
    for (const auto& chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (const auto& res : chain.residues) {
            if (res.seqid.num.value == seqnum) {
                return &res;
            }
        }
    }
    return nullptr;
}

// Helper: find atom by trimmed name in a gemmi residue
static const gemmi::Atom* find_atom(const gemmi::Residue& res, const std::string& atom_name_trimmed) {
    for (const auto& at : res.atoms) {
        std::string aname = at.name;
        size_t start = aname.find_first_not_of(" ");
        size_t end = aname.find_last_not_of(" ");
        std::string trimmed;
        if (start != std::string::npos) {
            trimmed = aname.substr(start, end - start + 1);
        }
        if (trimmed == atom_name_trimmed) {
            return &at;
        }
    }
    return nullptr;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, get_stats_around_ca) {
    molecules_container_t mc;
    int imol = mc.read_pdb("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(imol, 0) << "PDB load failed";

    int imap = mc.read_mtz("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.mtz",
                           "FWT", "PHWT", "", false, false);
    ASSERT_GE(imap, 0) << "MTZ load failed";

    clipper::Xmap<float> xmap = mc.get_xmap(imap);

    // Load PDB with gemmi for residue/atom access
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    // Constants (matching original test: n_steps=1, grid_box_radius=1.5f)
    const int n_steps = 1;
    const float grid_box_radius = 1.5f;

    // ---- Case 1: Standard residue with CB (//A/50, ARG) ----
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Residue* res = find_residue(model, "A", 50);
        ASSERT_NE(res, nullptr) << "residue not found";

        EXPECT_EQ(res->name, "ARG");

        const gemmi::Atom* ca_at = find_atom(*res, "CA");
        const gemmi::Atom* cb_at = find_atom(*res, "CB");
        const gemmi::Atom* c_at  = find_atom(*res, "C");
        const gemmi::Atom* n_at  = find_atom(*res, "N");

        ASSERT_NE(ca_at, nullptr);
        ASSERT_NE(c_at, nullptr);
        ASSERT_NE(n_at, nullptr);

        clipper::Coord_orth pt_ca = coot::co_gemmi(ca_at);
        clipper::Coord_orth pt_cb(0, 0, 0);
        if (cb_at) {
            pt_cb = coot::co_gemmi(cb_at);
        } else {
            pt_cb = clipper::Coord_orth(pt_ca.x() + 1.54, pt_ca.y(), pt_ca.z());
        }
        clipper::Coord_orth pt_c = coot::co_gemmi(c_at);
        clipper::Coord_orth pt_n = coot::co_gemmi(n_at);

        EXPECT_NEAR(pt_ca.x(), 8.709, 1e-3);
        EXPECT_NEAR(pt_ca.y(), -3.609, 1e-3);
        EXPECT_NEAR(pt_ca.z(), 56.949, 1e-3);

        auto axes = coot::make_axes_gemmi(pt_ca, pt_cb, pt_c, pt_n);
        ASSERT_FALSE(axes.empty()) << "axes construction failed";
        EXPECT_EQ(axes.size(), 3u);

        float step_size = grid_box_radius / static_cast<float>(n_steps);
        EXPECT_EQ(n_steps, 1);
        EXPECT_NEAR(step_size, 1.5f, 1e-4);

        auto result = coot::get_stats_around_ca_gemmi(res, axes, step_size, xmap, n_steps);
        double mean = std::get<0>(result);
        double var = std::get<1>(result);
        double mean_pos = std::get<2>(result);

        EXPECT_NEAR(mean, 0.274024, 1e-4);
        EXPECT_NEAR(var, 0.0328267, 1e-4);
        EXPECT_NEAR(mean_pos, 0.33681, 1e-4);
    }

    // ---- Case 2: Non-existent residue (//A/99999) ----
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Residue* res = find_residue(model, "A", 99999);
        EXPECT_EQ(res, nullptr) << "expected residue not found";
    }

    // ---- Case 3: Different residue on chain A (//A/100, PRO) ----
    {
        const gemmi::Model& model = st.models[0];
        const gemmi::Residue* res = find_residue(model, "A", 100);
        ASSERT_NE(res, nullptr) << "residue not found";

        EXPECT_EQ(res->name, "PRO");

        const gemmi::Atom* ca_at = find_atom(*res, "CA");
        const gemmi::Atom* cb_at = find_atom(*res, "CB");
        const gemmi::Atom* c_at  = find_atom(*res, "C");
        const gemmi::Atom* n_at  = find_atom(*res, "N");

        ASSERT_NE(ca_at, nullptr);
        ASSERT_NE(c_at, nullptr);
        ASSERT_NE(n_at, nullptr);

        clipper::Coord_orth pt_ca = coot::co_gemmi(ca_at);
        clipper::Coord_orth pt_cb(0, 0, 0);
        if (cb_at) {
            pt_cb = coot::co_gemmi(cb_at);
        } else {
            // pseudo CB for Glycine: place ~1.54A from CA
            pt_cb = clipper::Coord_orth(pt_ca.x() + 1.54, pt_ca.y(), pt_ca.z());
        }
        clipper::Coord_orth pt_c = coot::co_gemmi(c_at);
        clipper::Coord_orth pt_n = coot::co_gemmi(n_at);

        EXPECT_NEAR(pt_ca.x(), 48.663, 1e-3);
        EXPECT_NEAR(pt_ca.y(), 4.527, 1e-3);
        EXPECT_NEAR(pt_ca.z(), 53.538, 1e-3);

        auto axes = coot::make_axes_gemmi(pt_ca, pt_cb, pt_c, pt_n);
        ASSERT_FALSE(axes.empty()) << "axes construction failed";
        EXPECT_EQ(axes.size(), 3u);

        float step_size = grid_box_radius / static_cast<float>(n_steps);
        EXPECT_EQ(n_steps, 1);
        EXPECT_NEAR(step_size, 1.5f, 1e-4);

        auto result = coot::get_stats_around_ca_gemmi(res, axes, step_size, xmap, n_steps);
        double mean = std::get<0>(result);
        double var = std::get<1>(result);
        double mean_pos = std::get<2>(result);

        EXPECT_NEAR(mean, 0.391308, 1e-4);
        EXPECT_NEAR(var, 0.134556, 1e-4);
        EXPECT_NEAR(mean_pos, 0.528301, 1e-4);
    }
}