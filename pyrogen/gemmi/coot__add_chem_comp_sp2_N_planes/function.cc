#include "function.hh"
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/RDGeneral/Invariant.h>
#include <string>
#include <vector>
#include <cstdio>

using namespace RDKit;

void coot::add_chem_comp_sp2_N_planes_gemmi(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints) {

  typedef std::pair<std::string, double> d_pat;
  std::vector<d_pat> patterns;
  patterns.push_back(d_pat("[c,C][N^2;H2]([H])[H]", 0.02));  // N6 on Adenosine.
                                                             // Should the Hs be replaced by *s?
  int n_planes = 1; // counter for output text
  for (unsigned int ipat = 0; ipat < patterns.size(); ipat++) {
    RDKit::ROMol *query = SmartsToMol(patterns[ipat].first);
    std::vector<RDKit::MatchVectType> matches;
    bool recursionPossible = true;
    bool useChirality = true;
    bool uniquify = true;
    int matched = SubstructMatch(mol, *query, matches, uniquify, recursionPossible, useChirality);
    // int matched = false; // 20210923-PE FIXME
    if (false) // debug
      std::cout << "Matched " << matched << " sp2 N planes" << std::endl;
    for (unsigned int imatch = 0; imatch < matches.size(); imatch++) {
      if (matches[imatch].size() > 0) {
        std::cout << "matched sp2 N plane pattern: " << patterns[ipat].first << std::endl;
        std::string plane_id = "plane-sp2-N-";
        char s[100];
        snprintf(s, 99, "%d", n_planes);
        plane_id += std::string(s);
        try {
          std::vector<std::string> atom_names;
          for (unsigned int ii = 0; ii < matches[imatch].size(); ii++) {
            const RDKit::Atom *at_p = mol[matches[imatch][ii].second];

            // Unlike aromatics, the atoms of this type of plane
            // can be in more than one plane.

            std::string name = "";
            at_p->getProp("name", name);
            at_p->setProp("plane_id", plane_id);
            atom_names.push_back(name);
          }
          if (atom_names.size() > 3) {
            double dist_esd = patterns[ipat].second;
            dict_plane_restraint_t res(plane_id, atom_names, dist_esd);
            restraints->plane_restraint.push_back(res);
          }
        }
        catch (const KeyErrorException &kee) {

        }
        n_planes++;
      }
    }
    delete query;
  }
}