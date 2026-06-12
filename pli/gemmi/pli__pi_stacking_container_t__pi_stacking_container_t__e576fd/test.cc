#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include "function.hh"
#include "geometry/protein-geometry.hh"

// Helper: find residue by chain ID and sequence number
static gemmi::Residue* find_residue(const gemmi::Model& model, const std::string& chain_id, int seq_num) {
    for (const gemmi::Chain& chain : model.chains) {
        if (chain.name != chain_id) continue;
        for (const gemmi::Residue& res : chain.residues) {
            if (res.seqid.num.value == seq_num) {
                return const_cast<gemmi::Residue*>(&res);
            }
        }
    }
    return nullptr;
}

// Helper: collect all residues near a reference residue (simple distance check)
static std::vector<gemmi::Residue*> residues_near(
    gemmi::Residue* ref_res,
    const gemmi::Model& model,
    double radius)
{
    std::vector<gemmi::Residue*> result;
    if (!ref_res) return result;

    // Find the chain this residue belongs to
    const gemmi::Chain* ref_chain = nullptr;
    for (const gemmi::Chain& ch : model.chains) {
        for (const gemmi::Residue& r : ch.residues) {
            if (&r == ref_res) {
                ref_chain = &ch;
                break;
            }
        }
        if (ref_chain) break;
    }
    if (!ref_chain) return result;

    // Pick the first atom of the reference residue
    if (ref_res->atoms.empty()) return result;
    const gemmi::Vec3 ref_pos = ref_res->atoms[0].pos;

    double radius_sq = radius * radius;

    for (const gemmi::Chain& chain : model.chains) {
        for (const gemmi::Residue& res : chain.residues) {
            for (const gemmi::Atom& atom : res.atoms) {
                double dx = atom.pos.x - ref_pos.x;
                double dy = atom.pos.y - ref_pos.y;
                double dz = atom.pos.z - ref_pos.z;
                double dist_sq = dx*dx + dy*dy + dz*dz;
                if (dist_sq < radius_sq) {
                    result.push_back(const_cast<gemmi::Residue*>(&res));
                    break; // Add this residue once
                }
            }
        }
    }

    return result;
}

TEST(PiStackingGemmiTest, basic) {
    // Read structure
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_EQ(static_cast<int>(st.models.size()), 1);

    gemmi::Model& model = st.models[0];

    // Get protein geometry for monomer restraints
    coot::protein_geometry geom;
    geom.init_standard();

    // Case 1: Valid aromatic residue — PHE 80
    {
        gemmi::Residue* ref = find_residue(model, "A", 80);
        ASSERT_TRUE(ref != nullptr);
        EXPECT_EQ(ref->name, "PHE");

        // Get monomer restraints
        std::pair<bool, coot::dictionary_residue_restraints_t> p =
            geom.get_monomer_restraints_at_least_minimal("PHE", 0);
        ASSERT_TRUE(p.first);

        // Find nearby residues
        double residues_near_radius = 10.0;
        std::vector<gemmi::Residue*> residues = residues_near(ref, model, residues_near_radius);

        EXPECT_GE(residues.size(), 1u); // At least the reference residue

        // Empty aromatic ring list — mimics original constructor calling get_aromatic_ring_list
        // which returns a list specific to the residue type; with no matching rings, result is 0
        std::vector<std::vector<std::string>> aromatic_ring_list;

        // Construct pi_stacking_container via gemmi port
        pli::pi_stacking_container_t pi_stack =
            pli::pi_stacking_container_t_gemmi(p.second, residues, ref, aromatic_ring_list);

        EXPECT_EQ(pi_stack.size(), 0u);
    }

    // Case 2: Edge case — empty residues vector
    {
        gemmi::Residue* ref = find_residue(model, "A", 80);
        ASSERT_TRUE(ref != nullptr);

        std::pair<bool, coot::dictionary_residue_restraints_t> p =
            geom.get_monomer_restraints_at_least_minimal("PHE", 0);
        ASSERT_TRUE(p.first);

        std::vector<gemmi::Residue*> empty_residues;
        std::vector<std::vector<std::string>> aromatic_ring_list;
        pli::pi_stacking_container_t pi_stack_empty =
            pli::pi_stacking_container_t_gemmi(p.second, empty_residues, ref, aromatic_ring_list);
        EXPECT_EQ(pi_stack_empty.size(), 0u);
    }

    // Case 3: Ligand without monomer restraints (LZA 1299)
    {
        gemmi::Residue* ref = find_residue(model, "A", 1299);
        if (ref) {
            EXPECT_EQ(ref->name, "LZA");
            std::pair<bool, coot::dictionary_residue_restraints_t> p =
                geom.get_monomer_restraints_at_least_minimal("LZA", 0);
            EXPECT_FALSE(p.first);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
