#include "function.hh"

using namespace RDKit;

namespace coot {

void add_chem_comp_deloc_planes_gemmi(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints) {

   typedef std::pair<std::string, double> d_pat;

   std::vector<d_pat> patterns;
   patterns.push_back(d_pat("*C(=O)[O;H]",                 0.02));  // ASP carboxylate, valence model (H)
   patterns.push_back(d_pat("*[C;X3;^2](~O)~[O;X1]",       0.02));  // ASP carboxylate, no H
   patterns.push_back(d_pat("AC(=O)[N^2;H2,H1]([H])[A,H]", 0.02));  // ASN
   patterns.push_back(d_pat("*C(=N)[N^2;H2]([H])[A,H]",    0.02));  // amidine
   patterns.push_back(d_pat("CNC(=[NH])N([H])[H]",         0.02));  // guanidinium with H - testing
   patterns.push_back(d_pat("CNC(=[NH])N",                 0.02));  // guanidinium sans Hs
   patterns.push_back(d_pat("*[C;X3;^2](=O)[N;X3;^2;H1]([H])*", 0.02));  // amino

   // Martin's pattern, these should be weaker (than standard 0.02) though, I think
   patterns.push_back(d_pat("[*^2]=[*^2]-[*^2]=[*;X1;^2]", 0.04));

   int n_planes = 1;
   for (unsigned int ipat=0; ipat<patterns.size(); ipat++) {
      RDKit::ROMol *query = RDKit::SmartsToMol(patterns[ipat].first);
      std::vector<RDKit::MatchVectType>  matches;
      bool recursionPossible = true;
      bool useChirality = true;
      bool uniquify = true;
      int matched = RDKit::SubstructMatch(mol,*query,matches,uniquify,recursionPossible, useChirality);
      for (unsigned int imatch=0; imatch<matches.size(); imatch++) {
	 if (matches[imatch].size() > 0) {

	    if (true) { // debug
	       std::cout << "INFO:: matched deloc plane: " << patterns[ipat].first << " ";
	       std::cout << " (";
	       for (unsigned int iat=0; iat<matches[imatch].size(); iat++) {
		  unsigned int atom_idx = matches[imatch][iat].second;
		  try {
		     const RDKit::Atom *at_p = mol[atom_idx];
		     std::string atom_name;
		     at_p->getProp("name", atom_name);
		     std::cout << " " << atom_name;
		  }
		  catch (const KeyErrorException &kee) {
		     std::cout << " " << atom_idx;
		  }
	       }
	       std::cout << " )";
	       std::cout << std::endl;
	    }

	    std::string plane_id = "plane-deloc-";
	    char s[100];
	    snprintf(s,99,"%d", n_planes);
	    plane_id += std::string(s);
	    try {
	       std::vector<std::string> atom_names;
	       for (unsigned int ii=0; ii<matches[imatch].size(); ii++) {
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
		  coot::dict_plane_restraint_t res(plane_id, atom_names, dist_esd);
		  restraints->plane_restraint.push_back(res);
	       }
	    }

	    catch (const KeyErrorException &kee) {
	       std::cout << "ERROR:: in add_chem_comp_planes_deloc failed to get atom name"
			 << std::endl;
	    }
	    n_planes++;
	 }
      }
   }
}

}