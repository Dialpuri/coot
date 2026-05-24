#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/math.hpp>
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__omega_torsion/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__co/gemmi/function.hh"
#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__util__cis_trans_convert/gemmi/function.hh"
#include "function.hh"

namespace coot {
namespace util {

int
cis_trans_conversion_gemmi(gemmi::Residue *res_first, gemmi::Residue *res_second,
                           gemmi::Structure &mol, gemmi::Structure &standard_residues_mol) {
    int count = 0;

    // Build vectors of TNS (trans) and CIS (cis) residues from standard_residues_mol
    std::vector<gemmi::Residue*> trans_residues;
    std::vector<gemmi::Residue*> cis_residues;

    for (auto& model : standard_residues_mol.models) {
        for (auto& chain : model.chains) {
            for (auto& residue : chain.residues) {
                if (residue.name == "TNS") {
                    trans_residues.push_back(&residue);
                } else if (residue.name == "CIS") {
                    cis_residues.push_back(&residue);
                }
            }
        }
    }

    if (trans_residues.size() >= 2 && cis_residues.size() >= 2) {
        auto mol_residues_pair = std::make_pair(res_first, res_second);
        int istatus = coot::util::cis_trans_convert_gemmi(mol_residues_pair, trans_residues, cis_residues);
        if (istatus) count = 1;
    }

    return count;
}

} // namespace util
} // namespace coot