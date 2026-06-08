#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/elem.hpp>
#include "function.hh"

#include <utility>
#include <vector>
#include <string>

static void populate_restraints_from_gemmi(
    coot::restraints_container_t& rc,
    const gemmi::Model& model,
    const std::vector<std::pair<bool, const gemmi::Residue*>>& moving_residues)
{
    rc.residues_vec_moving_set.clear();
    rc.from_residue_vector = !moving_residues.empty();
    for (const auto& pr : moving_residues) {
        if (pr.second) {
            coot::restraints_container_t::residue_identifier rid;
            rid.seqnum = pr.second->seqid.num.value;
            rid.icode  = pr.second->seqid.icode;
            rc.residues_vec_moving_set.insert(rid);
        }
    }

    rc.atom.clear();
    rc.residue.clear();
    rc.chain_ids.clear();

    for (const auto& pr : moving_residues) {
        if (!pr.second) continue;
        const gemmi::Residue* res = pr.second;
        std::string chain_name;
        for (const auto& chain : model.chains) {
            for (const auto& r : chain.residues) {
                if (&r == res) {
                    chain_name = chain.name;
                    break;
                }
            }
            if (!chain_name.empty()) break;
        }
        for (const auto& atom : res->atoms) {
            rc.atom.push_back(const_cast<gemmi::Atom*>(&atom));
            rc.residue.push_back(const_cast<gemmi::Residue*>(res));
            rc.chain_ids.push_back(chain_name);
        }
    }
    rc.n_atoms = rc.atom.size();
    rc.istart_res = 0;
    rc.iend_res   = 1000000;
}

static void populate_all_atoms(
    coot::restraints_container_t& rc,
    const gemmi::Model& model)
{
    rc.atom.clear();
    rc.residue.clear();
    rc.chain_ids.clear();

    for (const auto& chain : model.chains) {
        for (const auto& res : chain.residues) {
            for (const auto& atom : res.atoms) {
                rc.atom.push_back(const_cast<gemmi::Atom*>(&atom));
                rc.residue.push_back(const_cast<gemmi::Residue*>(&res));
                rc.chain_ids.push_back(chain.name);
            }
        }
    }
    rc.n_atoms = rc.atom.size();
    rc.istart_res = 0;
    rc.iend_res   = 1000000;
}

TEST(RestraintsContainerTest, init_shared_post_gemmi_basic) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    ASSERT_FALSE(st.models.empty());
    const gemmi::Model& model = st.models[0];

    const gemmi::Residue* first_res = nullptr;
    for (const auto& chain : model.chains) {
        if (chain.name == "A" && !chain.residues.empty()) {
            first_res = &chain.residues[0];
            break;
        }
    }
    ASSERT_NE(first_res, nullptr);

    EXPECT_EQ(first_res->seqid.num.value, 1);
    EXPECT_EQ(first_res->name, "MET");
    EXPECT_EQ(first_res->atoms.size(), 19u);

    std::vector<std::pair<bool, const gemmi::Residue*>> residues;
    residues.push_back(std::make_pair(false, first_res));

    coot::restraints_container_t rc;
    rc.do_hydrogen_atom_refinement = false;
    rc.cryo_em_mode = false;
    rc.do_neutron_refinement = false;

    populate_restraints_from_gemmi(rc, model, residues);

    EXPECT_EQ(rc.n_atoms, 19u);
    EXPECT_EQ(rc.from_residue_vector, true);

    rc.init_shared_post_gemmi({});

    EXPECT_EQ(rc.initial_position_params_vec.size(), static_cast<size_t>(3 * rc.n_atoms));
    EXPECT_EQ(rc.bonded_atom_indices.size(), rc.n_atoms);
    EXPECT_EQ(rc.use_map_gradient_for_atom.size(), rc.n_atoms);

    EXPECT_EQ(rc.udd_bond_angle, -1);
    EXPECT_EQ(rc.udd_atom_index_handle, -1);

    EXPECT_GE(rc.initial_position_params_vec.size(), 3u);
    std::cout << "Atom 0 pos: "
              << rc.initial_position_params_vec[0] << " "
              << rc.initial_position_params_vec[1] << " "
              << rc.initial_position_params_vec[2] << std::endl;

    EXPECT_NEAR(rc.initial_position_params_vec[0], 12.334, 1e-2);
    EXPECT_NEAR(rc.initial_position_params_vec[1], 2.772, 1e-2);
    EXPECT_NEAR(rc.initial_position_params_vec[2], 74.951, 1e-2);
}

TEST(RestraintsContainerTest, init_shared_post_gemmi_all_chain_a) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    ASSERT_FALSE(st.models.empty());
    const gemmi::Model& model = st.models[0];

    std::vector<std::pair<bool, const gemmi::Residue*>> moving_residues;
    for (const auto& chain : model.chains) {
        if (chain.name != "A") continue;
        for (const auto& res : chain.residues) {
            moving_residues.push_back(std::make_pair(false, &res));
        }
    }
    ASSERT_GT(moving_residues.size(), 0u);

    coot::restraints_container_t rc;
    rc.do_hydrogen_atom_refinement = false;
    rc.cryo_em_mode = false;
    rc.do_neutron_refinement = false;

    populate_restraints_from_gemmi(rc, model, moving_residues);

    EXPECT_GT(rc.n_atoms, 0u);
    EXPECT_EQ(rc.from_residue_vector, true);

    rc.init_shared_post_gemmi({});

    EXPECT_EQ(rc.initial_position_params_vec.size(), static_cast<size_t>(3 * rc.n_atoms));
    EXPECT_EQ(rc.bonded_atom_indices.size(), rc.n_atoms);
    EXPECT_EQ(rc.use_map_gradient_for_atom.size(), rc.n_atoms);

    int non_h_with_grad = 0;
    for (size_t i = 0; i < rc.n_atoms; i++) {
        bool is_h = rc.atom[i]->element.is_hydrogen();
        if (!is_h && rc.use_map_gradient_for_atom[i]) non_h_with_grad++;
    }
    EXPECT_GT(non_h_with_grad, 0);
}

TEST(RestraintsContainerTest, init_shared_post_gemmi_with_fixed_atoms) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    ASSERT_FALSE(st.models.empty());
    const gemmi::Model& model = st.models[0];

    const gemmi::Residue* first_res = nullptr;
    for (const auto& chain : model.chains) {
        if (chain.name == "A" && !chain.residues.empty()) {
            first_res = &chain.residues[0];
            break;
        }
    }
    ASSERT_NE(first_res, nullptr);

    std::vector<std::pair<bool, const gemmi::Residue*>> residues;
    residues.push_back(std::make_pair(false, first_res));

    coot::restraints_container_t rc;
    rc.do_hydrogen_atom_refinement = false;
    rc.cryo_em_mode = false;
    rc.do_neutron_refinement = false;

    populate_restraints_from_gemmi(rc, model, residues);

    std::vector<coot::atom_spec_t> fixed_specs;
    std::string chain_name = "A";
    int count = 0;
    for (const auto& atom : first_res->atoms) {
        if (count >= 3) break;
        coot::atom_spec_t spec(atom, *first_res, chain_name);
        fixed_specs.push_back(spec);
        count++;
    }

    EXPECT_FALSE(fixed_specs.empty());

    rc.init_shared_post_gemmi(fixed_specs);

    EXPECT_GT(rc.fixed_atom_indices.size(), 0u);

    for (int idx : rc.fixed_atom_indices) {
        EXPECT_FALSE(rc.use_map_gradient_for_atom[idx]);
    }
}

TEST(RestraintsContainerTest, init_shared_post_gemmi_non_residue_vector) {
    gemmi::Structure st = gemmi::read_pdb_file(
        "/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    ASSERT_FALSE(st.models.empty());
    const gemmi::Model& model = st.models[0];

    coot::restraints_container_t rc;
    rc.do_hydrogen_atom_refinement = false;
    rc.cryo_em_mode = false;
    rc.do_neutron_refinement = false;
    rc.from_residue_vector = false;

    populate_all_atoms(rc, model);

    EXPECT_GT(rc.n_atoms, 0u);

    rc.init_shared_post_gemmi({});

    EXPECT_EQ(rc.initial_position_params_vec.size(), static_cast<size_t>(3 * rc.n_atoms));
    EXPECT_EQ(rc.bonded_atom_indices.size(), rc.n_atoms);
    EXPECT_EQ(rc.use_map_gradient_for_atom.size(), rc.n_atoms);

    int non_h_with_grad = 0;
    for (size_t i = 0; i < rc.n_atoms; i++) {
        if (!rc.atom[i]->element.is_hydrogen() && rc.use_map_gradient_for_atom[i]) {
            non_h_with_grad++;
        }
    }
    EXPECT_GT(non_h_with_grad, 0);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}