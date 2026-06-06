#pragma once

#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__cis_trans_convert/gemmi/function.hh"

namespace coot {
namespace util {

int cis_trans_conversion_gemmi(
    gemmi::CRA cra,
    bool is_N_flag,
    gemmi::Model* mol,
    gemmi::Model* standard_residues_model)
{
    if (!cra.atom || !cra.residue || !cra.chain) return 0;
    if (!mol || !standard_residues_model) return 0;

    int offset = is_N_flag ? -1 : 0;
    int resno_1 = cra.residue->seqid.num.value + offset;
    int resno_2 = resno_1 + 1;
    std::string chain_id = cra.chain->name;

    // Find the two consecutive residues in mol
    std::vector<gemmi::Residue*> mol_residues;
    for (auto& chain : mol->chains) {
        if (chain.name != chain_id) continue;
        for (auto& res : chain.residues) {
            int seq = res.seqid.num.value;
            if (seq == resno_1 || seq == resno_2) {
                mol_residues.push_back(&res);
            }
        }
    }

    if (mol_residues.size() < 2) {
        std::cerr << "ERROR: failed to get mol residues in cis_trans_convert" << std::endl;
        return 0;
    }

    // Sort by sequence number to ensure correct order
    std::sort(mol_residues.begin(), mol_residues.end(), [](gemmi::Residue* a, gemmi::Residue* b) {
        return a->seqid.num.value < b->seqid.num.value;
    });

    // Find TNS residues from standard residues model
    std::vector<gemmi::Residue*> trans_residues;
    for (auto& chain : standard_residues_model->chains) {
        for (auto& res : chain.residues) {
            if (res.name == "TNS") {
                trans_residues.push_back(&res);
            }
        }
    }

    if (trans_residues.size() < 2) {
        std::cerr << "ERROR: failed to get trans residues in cis_trans_convert" << std::endl;
        return 0;
    }

    // Find CIS residues from standard residues model
    std::vector<gemmi::Residue*> cis_residues;
    for (auto& chain : standard_residues_model->chains) {
        for (auto& res : chain.residues) {
            if (res.name == "CIS") {
                cis_residues.push_back(&res);
            }
        }
    }

    if (cis_residues.size() < 2) {
        std::cerr << "ERROR: failed to get cis residues in cis_trans_convert " << cis_residues.size() << std::endl;
        return 0;
    }

    std::pair<gemmi::Residue*, gemmi::Residue*> mol_residues_pair(mol_residues[0], mol_residues[1]);
    return cis_trans_convert_gemmi(mol_residues_pair, trans_residues, cis_residues);
}

} // namespace util
} // namespace coot