#include "function.hh"

#include <iostream>
#include <cmath>
#include <vector>

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

#include <gemmi/model.hpp>

namespace {

// Build a clipper::Atom_list from a gemmi::Model
clipper::Atom_list atoms_from_gemmi_model(gemmi::Model *model) {
    std::vector<clipper::Atom> atoms;
    // Count atoms first
    size_t n = 0;
    for (gemmi::Chain &chain : model->chains)
        for (gemmi::Residue &res : chain.residues)
            n += res.atoms.size();
    atoms.reserve(n);

    for (gemmi::Chain &chain : model->chains) {
        for (gemmi::Residue &res : chain.residues) {
            for (gemmi::Atom &atom : res.atoms) {
                clipper::Atom clipper_atom;
                // Element from gemmi — Element.name() returns "C", "O", etc.
                clipper_atom.set_element(atom.element.name());
                // Position (clipper uses Coord_orth which is x,y,z)
                clipper_atom.set_coord_orth(
                    clipper::Coord_orth(atom.pos.x, atom.pos.y, atom.pos.z));
                // Occupancy and B-factor (u_iso = B / 8pi^2)
                clipper_atom.set_occupancy(atom.occ);
                clipper_atom.set_u_iso(atom.b_iso / (8.0 * M_PI * M_PI));

                atoms.push_back(std::move(clipper_atom));
            }
        }
    }

    return clipper::Atom_list(atoms);
}

} // anonymous namespace

void coot::util::sfcalc_genmap_gemmi(
    gemmi::Model *mol,
    const clipper::HKL_data<clipper::data32::F_sigF> &fobs_in,
    const clipper::HKL_data<clipper::data32::Flag> &free,
    clipper::Xmap<float> *xmap_p) {

   if (fobs_in.num_obs() == 0) {
      std::cout << "sfcalc_genmap_gemmi(): No Fobs reflections\n";
      return;
   }
   enum ANISO { NONE, FOBS, FCAL };
   bool bulk  = true;
   ANISO aniso = NONE;
   int freeflag = 0;
   int n_refln = 1000;
   int n_param = 20;
   int verbose = 0;
   typedef clipper::HKL_data_base::HKL_reference_index HRI;

   // modifyable fobs (because we want to anisotropically scale Fobs)
   clipper::HKL_data<clipper::data32::F_sigF> fobs(fobs_in);

   // get a list of all the atoms from gemmi model
   clipper::Atom_list atoms = atoms_from_gemmi_model(mol);

   clipper::HKL_info hkls;
   hkls.init(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling(), true);
   double bulkfrc, bulkscl;

   // calculate structure factors
   clipper::HKL_data<clipper::data32::F_phi> fc(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling());
   if (bulk) {
      clipper::SFcalc_obs_bulk<float> sfcb;
      sfcb(fc, fobs, atoms);
      bulkfrc = sfcb.bulk_frac();
      bulkscl = sfcb.bulk_scale();
   } else {
      clipper::SFcalc_aniso_fft<float> sfc;
      sfc(fc, atoms);
      bulkfrc = bulkscl = 0.0;
   }

  // do anisotropic scaling
   if (aniso != NONE)  {
      clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
      clipper::SFscale_aniso<float> sfscl;
      if ( aniso == FOBS ) sfscl(fobs, fc);
      if ( aniso == FCAL ) sfscl(fc, fobs);
      std::cout << "\nAnisotropic scaling:\n" << sfscl.u_aniso_orth(F).format() << "\n";
   }

   // now do sigmaa calc
   clipper::HKL_data<clipper::data32::F_phi> f_best(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling());
   clipper::HKL_data<clipper::data32::F_phi> f_diff(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling());
   clipper::HKL_data<clipper::data32::Phi_fom> phiw(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling());
   clipper::HKL_data<clipper::data32::Flag>    flag(fobs.spacegroup(),fobs.cell(),fobs.hkl_sampling());

   for ( HRI ih = flag.first(); !ih.last(); ih.next() )
      if (!fobs[ih].missing() && (free[ih].missing()||free[ih].flag()==freeflag))
         flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
      else
         flag[ih].flag() = clipper::SFweight_spline<float>::NONE;

   // do sigmaa calc
   clipper::SFweight_spline<float> sfw(n_refln, n_param);
   sfw(f_best, f_diff, phiw, fobs, fc, flag);

   if (true) {
      unsigned int n_nans_fobs   = 0;
      unsigned int n_nans_fc     = 0;
      unsigned int n_nans_f_diff = 0;
      for (HRI ih = fobs.first(); !ih.last(); ih.next()) {
         if (!fobs[ih].missing() && (free[ih].missing()||free[ih].flag()==freeflag)) {
            if (clipper::Util::isnan(fobs[ih].f())) n_nans_fobs++;
            if (clipper::Util::isnan(fc[ih].f()))   n_nans_fc++;
            if (clipper::Util::isnan(f_diff[ih].f())) n_nans_f_diff++;
         }
      }
      std::cout << "DEBUG:: sfcalc_genmap() the nan count: " << n_nans_fobs << " " << n_nans_fc << " " << n_nans_f_diff << std::endl;
   }

   // calc abcd (needed?)
   clipper::HKL_data<clipper::data32::ABCD> abcd(hkls);
   abcd.compute(phiw, clipper::data32::Compute_abcd_from_phifom());

   // now calc R and R-free
   std::vector<double> params(n_param, 1.0);
   clipper::BasisFn_spline basisfn(fobs, n_param, 1.0);
   clipper::TargetFn_scaleF1F2<clipper::data32::F_phi,clipper::data32::F_sigF> targetfn(fc, fobs);
   clipper::ResolutionFn rfn(hkls, basisfn, targetfn, params);
   double r1w, f1w, r1f, f1f, Fo, Fc;
   r1w = f1w = r1f = f1f = 0.0;
   for (HRI ih = fobs.first(); !ih.last(); ih.next()) {
      if (!fobs[ih].missing()) {
         Fo = fobs[ih].f();
         Fc = sqrt( rfn.f(ih) ) * fc[ih].f();
         if (free[ih].flag() == freeflag) {
            r1f += fabs(Fo - Fc);
            f1f += Fo;
         } else {
            r1w += fabs(Fo - Fc);
            f1w += Fo;
         }
      }
   }
   r1f /= clipper::Util::max(f1f, 0.1);
   r1w /= clipper::Util::max(f1w, 0.1);
   if (false)
      std::cout << "\n R-factor      : " << r1w << "\n Free R-factor : " << r1f << "\n";

   // now make a map
   xmap_p->fft_from(f_diff);

   // DIAGNOSTIC OUTPUT
   if (false) {
      std::cout << "\n Bulk Correction Volume: " << bulkfrc;
      std::cout << "\n Bulk Correction Factor: " << bulkscl << "\n";
      std::cout << "\nNumber of spline params: " << sfw.params_scale().size() << "\n";
      clipper::BasisFn_spline basisfn2( hkls, sfw.params_scale().size(), 1.0 );
      printf("\n $TABLE: Sigmaa statistics :\n $GRAPHS:scale vs resolution:N:1,2:\n        :lack of closure vs resolution:N:1,3:\n $$\n 1/resol^2   scale   lack_of_closure $$\n $$\n");
      for (int i = 0; i <= 20; i++) {
         double s = hkls.resolution().invresolsq_limit()*double(i)/20.0;
         printf("%6.3f %12.3f %12.3f\n", s, basisfn2.f_s(s,sfw.params_scale()),
                                         basisfn2.f_s(s,sfw.params_error()));
      }
      printf(" $$\n");
   }
}