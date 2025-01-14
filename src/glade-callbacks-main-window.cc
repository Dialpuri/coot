/* src/glade-callbacks.cc
 *
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * Copyright 2008 The University of Oxford
 * Copyright 2015, 2016 by Medical Research Council
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */


#include "Python.h"

#include <iostream>
#include <gtk/gtk.h>

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "coot-fileselections.h"
#include "positioned-widgets.h"
#include "interface.h"
#include "coot-references.h"

// put preferences functions into their own file, not here.
#include "coot-preferences.h"
#include "c-interface-preferences.h"
#include "rotate-translate-modes.hh"
#include "restraints-editor-c.h"
#include "generic-display-objects-c.h"
#include "c-interface-refmac.h"
#include "gtk-widget-conversion-utils.h"
#include "curlew.h"
#include "read-phs.h"
#include "gtk-manual.h"
#include "c-interface-refine.h"
#include "widget-from-builder.hh"
#include "read-molecule.hh" // 20230621-PE now with std::string args

// this from callbacks.h (which I don't want to include here)
typedef const char entry_char_type;

#include <vector>
#include "utils/coot-utils.hh"
#include "graphics-info.h"

#include "cc-interface.hh" // for read_ccp4_map()


// Let's put the new refinement and regularization control tools together here
// (although not strictly main window)


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_refine_control_button_clicked(GtkButton       *button,
                                               gpointer         user_data) {
   graphics_info_t::show_refinement_and_regularization_parameters_frame();
}




extern "C" G_MODULE_EXPORT
void
on_model_toolbar_select_map_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                           G_GNUC_UNUSED gpointer         user_data) {

   show_select_map_dialog();
}

// <signal name="toggled" handler="on_model_toolbar_range_define_togglebutton_toggled" swapped="no"/>

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_range_define_togglebutton_toggled(GtkToggleButton *togglebutton,
                                                   G_GNUC_UNUSED gpointer user_data) {

   gboolean active = gtk_toggle_button_get_active(togglebutton);
   graphics_info_t g;
   if (active) {
      if (g.in_range_define == 0) g.in_range_define = 1;
   } else {
      if (g.in_range_define == 1) g.in_range_define = 0;
      if (g.in_range_define == 2) g.in_range_define = 0;
   }
   std::cout << "here now with active " << active << " in_range_define " << g.in_range_define << std::endl;

}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_fixed_atoms_menubutton_activated(GtkMenuButton *button,
                                                  gpointer       user_data) {
   GtkWidget *w = wrapped_create_fixed_atom_dialog();
   gtk_widget_set_visible(w, TRUE);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rigid_body_fit_togglebutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                        gpointer         user_data) {

   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active) {
      printf("Rigid Body:\n");
      do_rigid_body_refine(1);
   } else {
      do_rigid_body_refine(0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rot_trans_togglebutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                        gpointer         user_data)
{
   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active) {
      do_rot_trans_setup(1);
   } else {
      do_rot_trans_setup(0);
   }
}

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_auto_fit_rotamer_button_clicked(GtkButton *button,
                                                 gpointer   user_data) {

#if 0 // 20220813-PE there is no setup now
   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active)
      setup_auto_fit_rotamer(1);
   else
      setup_auto_fit_rotamer(0);
#endif

   graphics_info_t g;
   std::pair<int, mmdb::Atom *> aa = g.get_active_atom();
   int imol = aa.first;
   if (is_valid_model_molecule(imol)) {
      std::string alt_conf = aa.second->altLoc;
      coot::residue_spec_t res_spec(coot::atom_spec_t(aa.second));
      g.auto_fit_rotamer_ng(imol, res_spec, alt_conf);
   }

}


// 20220812-PE this used to be a toggle-button - it expected an atom pick after toggling the
// button. These days we find the rotamers for the residue at the centre of the screen.
// There is no longer any "setup_rotamers()"
//
extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rotamers_button_clicked(GtkButton *button,
                                         gpointer  user_data) {

   // Find rotamers for the residue at the centre of the screen
   graphics_info_t g;
   std::pair<int, mmdb::Atom *> aa = g.get_active_atom();
   int imol = aa.first;
   if (is_valid_model_molecule(imol)) {
      g.do_rotamers(imol, aa.second);
   }
}

// delete this when done
extern "C" G_MODULE_EXPORT
void
on_model_toolbar_edit_chi_angles_togglebutton_toggled(GtkToggleButton *toggletoolbutton,
                                                                          gpointer         user_data) {

   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active) {
      setup_edit_chi_angles(1);
   } else {
      setup_edit_chi_angles(0);
      set_show_chi_angle_bond(0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_edit_chi_angles_button_clicked(GtkButton *button,
                                                gpointer   user_data) {

   graphics_info_t g;
   auto active_atom = g.get_active_atom();
   int imol = active_atom.first;
   if (is_valid_model_molecule(imol)) {
      auto &m = g.molecules[imol];
      mmdb:: Atom *at = active_atom.second;
      int idx = m.get_atom_index(at);
      g.execute_edit_chi_angles(idx, imol);
   }
}



extern "C" G_MODULE_EXPORT
void
on_model_toolbar_flip_peptide_button_clicked(GtkButton *button,
                                             gpointer   user_data) {
   graphics_info_t g;
   // we want the actual atom, not the CA of the residue
   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa = g.active_atom_spec_simple();
   if (aa.first) {
      int imol = aa.second.first;
      if (is_valid_model_molecule(imol)) {
         auto &m = g.molecules[imol];
         coot::atom_spec_t atom_spec(aa.second.second);
         m.pepflip(atom_spec);
         g.graphics_draw();
      }
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_side_chain_180_button_clicked(GtkButton *button,
                                               gpointer   user_data) {

#if 0 // do it directly these days
   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
  if (active)
    setup_180_degree_flip(1);
  else
    setup_180_degree_flip(0);
#endif

   // look at check_if_in_180_degree_flip_define, it checks for
   // intermediate atoms - we should do that here too.

   std::cout << "Here in on_model_toolbar_side_chain_180_button_clicked()" << std::endl;
   graphics_info_t g;
   auto active_atom = g.get_active_atom();
   int imol = active_atom.first;
   std::cout << "imol = " << imol << " atom = " << active_atom.second << std::endl;
   if (is_valid_model_molecule(imol)) {
      mmdb::Atom *at = active_atom.second;
      std::string alt_conf(at->altLoc);
      coot::atom_spec_t at_spec(at);
      coot::residue_spec_t spec(at_spec);
      auto &m = g.molecules[imol];
      // change this signature to use a residue spec and an alt_conf.
      int istatus = m.do_180_degree_side_chain_flip(spec.chain_id, spec.res_no, spec.ins_code, alt_conf, g.Geom_p());
      g.graphics_draw();
   }

}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_mutate_and_autofit_menubutton_activate(GtkMenuButton *menubutton,
                                                        gpointer       user_data) {

   // this function seems not to be called on menu button click. Hmmm.
   std::cout << "on_model_toolbar_mutate_and_autofit_menubutton_active "
             << " select the right menu for mutate_and_autofit_menubutton here" << std::endl;
}


extern "C" G_MODULE_EXPORT
void
on_simple_mutate_menubutton_activate(GtkMenuButton *menubutton,
                                     gpointer       user_data) {

   // "activate" happens when the popover menu is shown. That is not what we want.

   std::cout << "on_simple_mutate_menubutton_activate "
             << " select the right menu for mutate_and_autofit_menubutton here" << std::endl;
}

#include "setup-gui-components.hh"

// try again with just a button

extern "C" G_MODULE_EXPORT
void
on_simple_mutate_button_clicked(GtkButton *button,
                                gpointer   user_data) {

   GtkWidget *menu_button = widget_from_builder("simple_mutate_menubutton");
   if (menu_button) {
      graphics_info_t g;
      auto active_atom = g.get_active_atom();
      int imol = active_atom.first;
      if (is_valid_model_molecule(imol)) {
         gtk_widget_set_visible(menu_button, TRUE);
         mmdb::Atom *atom = active_atom.second;
         mmdb::Residue *residue_p = atom->residue;
         if (residue_p) {
            if (coot::util::is_nucleotide_by_dict(residue_p, *g.Geom_p())) {
               add_typed_menu_to_mutate_menubutton("SIMPLE", "NUCLEIC-ACID");
            } else {
               add_typed_menu_to_mutate_menubutton("SIMPLE", "PROTEIN");
            }
            g_signal_emit_by_name(menu_button, "activate", NULL);
         }
      }
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_add_terminal_residue_button_clicked(GtkButton *button,
                                                     gpointer   user_data)
{
   graphics_info_t g;
   auto active_atom = g.get_active_atom();
   int imol = active_atom.first;
   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = active_atom.second->residue;
      coot::atom_spec_t atom_spec(active_atom.second);
      std::string chain_id = atom_spec.chain_id;
      mmdb::Atom *atom = active_atom.second;
      std::string terminus_type = g.molecules[imol].get_term_type(atom);

      if (coot::util::is_nucleotide_by_dict_dynamic_add(residue_p, g.Geom_p())) {
         std::cout << "add nucleotide" << std::endl;
         g.execute_simple_nucleotide_addition(imol, terminus_type, residue_p, chain_id);
      } else {
         coot::residue_spec_t residue_spec(atom_spec);
         std::string new_type = "ALA";
         g.execute_add_terminal_residue(imol, terminus_type, residue_p, chain_id, new_type, true);
      }
   }

}




extern "C" G_MODULE_EXPORT
void
on_model_toolbar_add_alt_conf_button_clicked(GtkButton *button,
                                             gpointer   user_data) {

   // altconf(); 20240304-PE this was this old "dialog before atom pick" method

   graphics_info_t g;
   auto active_atom_pair = g.get_active_atom();
   int imol = active_atom_pair.first;
   if (is_valid_model_molecule(imol)) {
      // split_residue_range(imol, add_alt_conf_atom_index, naii.atom_index);  // from graphics-info-defines
      // For today, we will make this just a single residue.
      const auto &active_atom = active_atom_pair.second;
      coot::atom_spec_t aas(active_atom);
      g.split_residue(imol, aas.chain_id, aas.res_no, aas.ins_code, aas.alt_conf);
   }
}




extern "C" G_MODULE_EXPORT
void
on_model_toolbar_undo_button_clicked(GtkButton       *button,
                                     gpointer         user_data) {
   apply_undo();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_redo_button_clicked(GtkButton       *button,
                                     gpointer         user_data) {
  apply_redo();
}


extern "C" G_MODULE_EXPORT
void
on_save_coordinates_filechooser_dialog_response(GtkDialog       *dialog,
                                                gint             response_id,
                                                gpointer         user_data) {
   if (response_id == GTK_RESPONSE_OK) {

      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *fnc = g_file_info_get_name(file_info);
      if (fnc) {

         // imol set in
         // on_save_coords_dialog_save_button_clicked(GtkButton       *button,
         //                                           gpointer         user_data)
         int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "imol"));
         save_coordinates(imol, fnc);
      }
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }
}

extern "C" G_MODULE_EXPORT
void
on_dataset_filechooser_dialog_response(GtkDialog       *dialog,
                                       gint             response_id,
                                       gpointer         user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      // const char *fnc = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *fnc = g_file_info_get_name(file_info);
      if (fnc) {
         int is_auto = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "is_auto"));
         std::string fn(fnc);
         if (is_auto == 1) {
            auto_read_make_and_draw_maps_from_mtz(fnc);
         } else {
            manage_column_selector(fnc); // move the function declaration into a c++ header one day
         }
      }
      save_directory_from_filechooser(GTK_WIDGET(dialog));
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }
}

extern "C" G_MODULE_EXPORT
void
on_map_name_filechooser_dialog_response(GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      // const char *fnc = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *fnc = g_file_info_get_name(file_info);
      if (fnc) {
         std::string fn(fnc);
         // std::cout << "Now do something with " << fn << std::endl;
         bool is_diff_map = false;
         GtkWidget *checkbutton = widget_from_builder("map_filechooser_is_difference_map_button");
         // 20220809-PE GTK4 (post merge) - this is a checkbutton FIXME
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton)))
            is_diff_map = true;
         std::cout << "----------------- read_ccp4_map() " << std::endl;
         read_ccp4_map(fn, is_diff_map);
         std::cout << "----------------- done read_ccp4_map() " << std::endl;
      }
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }

   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   
}


#include "rsr-functions.hh"


extern "C" G_MODULE_EXPORT
void
on_add_an_atom_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_water_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                    G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Water");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_lithium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                      G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Li");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_sodium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Na");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_potasium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                       G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("K");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_magnesium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                        G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Mg");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_calcium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                      G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Ca");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_strontium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                        G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Sr");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_nickel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Ni");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_copper_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Cu");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_zinc_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                   G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Zn");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_chlorine_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                       G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Cl");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_bromine_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                      G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Br");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_iodine_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("I");
   gtk_widget_set_visible(box, FALSE);
}

 
extern "C" G_MODULE_EXPORT
void
on_get_monomer_ok_button_clicked(GtkButton       *button,
                                 gpointer         user_data) {

   GtkWidget *entry = widget_from_builder("get_monomer_entry");
   if (entry) {
      handle_get_monomer_code(entry);
   }
   GtkWidget *frame = widget_from_builder("get_monomer_frame");
   gtk_widget_set_visible(frame, FALSE);
}


extern "C" G_MODULE_EXPORT
void on_generic_overlay_frame_cancel_button_clicked(GtkButton       *button,
                                                    gpointer         user_data) {
   GtkWidget* frame_widget = GTK_WIDGET(user_data);
   if(frame_widget) {
      gtk_widget_set_visible(frame_widget, FALSE);
   } else {
      g_error("ERROR:: in on_generic_overlay_frame_cancel_button_clicked() 'user_data' is NULL. Cannot hide overlay frame.");
   }
}

// 20240518-PE this is the only function in this file that uses python. Hmm.
// So rewrite mutate_by_overlap() into C++ (non-trivial)
//
#include "cc-interface-scripting.hh"

extern "C" G_MODULE_EXPORT
void on_replace_residue_ok_button_clicked(GtkButton *button,
                                          gpointer user_data) {

   GtkWidget *frame = widget_from_builder("replace_residue_frame");
   GtkWidget *entry = widget_from_builder("replace_residue_entry");
   std::string new_residue_type = gtk_editable_get_text(GTK_EDITABLE(entry));

   graphics_info_t g;
   std::pair<int, mmdb::Atom *> aa = g.get_active_atom();
   int imol = aa.first;
   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = aa.second->residue;
      if (residue_p) {
         coot::protein_geometry &geom = *graphics_info_t::Geom_p();
         std::pair<bool, coot::dictionary_residue_restraints_t> rp = geom.get_monomer_restraints(new_residue_type, imol);
         if (rp.first) {
            const auto &restraints_new_type = rp.second;
            mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
            std::string current_residue_type = residue_p->GetResName();
            std::pair<bool, coot::dictionary_residue_restraints_t> rp_current =
               geom.get_monomer_restraints(current_residue_type, imol);
            if (rp_current.first) {
               const auto &restraints_current_type = rp_current.second;
               int status = coot::util::mutate_by_overlap(residue_p, mol, restraints_current_type, restraints_new_type);
               if (status == 0)
                  graphics_info_t::log.log(logging::WARNING, "mutate_by_overlap() failed");
            }
         }
      }
   }
   gtk_widget_set_visible(frame, FALSE);

}

extern "C" G_MODULE_EXPORT
void on_replace_residue_cancel_button_clicked(GtkButton *button,
                                              gpointer user_data) {

   GtkWidget *frame = widget_from_builder("replace_residue_frame");
   gtk_widget_set_visible(frame, FALSE);

}

#include "c-interface-ligands.hh"

extern "C" G_MODULE_EXPORT
void
on_smiles_to_simple_3d_ok_button_clicked(GtkButton       *button,
                                         gpointer         user_data) {

   GtkWidget *entry = widget_from_builder("smiles_to_simple_3d_entry");
   if (entry) {
      std::string smiles_string = gtk_editable_get_text(GTK_EDITABLE(entry));
      smiles_to_simple_3d(smiles_string);
   }
   GtkWidget *frame = widget_from_builder("smiles_to_simple_3d_frame");
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_diff_map_peaks_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *vbox = widget_from_builder("dialog-vbox78");
   clear_diff_map_peaks();
   gtk_widget_set_visible(vbox, FALSE);
   graphics_info_t::hide_vertical_validation_frame_if_appropriate();
}

extern "C" G_MODULE_EXPORT
void
on_diff_map_peaks_update_button_clicked(GtkButton *button,
                                       gpointer         user_data) {
   graphics_info_t g;
   g.fill_difference_map_peaks_button_box(); 

}

#include "dynamic-validation.hh"

extern "C" G_MODULE_EXPORT
void
on_dynamic_validation_update_button_clicked(GtkButton *button,
                                            gpointer  user_data) {

   update_dynamic_validation();
}

extern "C" G_MODULE_EXPORT
void
on_dynamic_validation_close_button_clicked(GtkButton *button,
                                           gpointer  user_data) {

   GtkWidget *validation_vbox = widget_from_builder("dynamic_validation_vbox");
   gtk_widget_set_visible(validation_vbox, FALSE);
   graphics_info_t::hide_vertical_validation_frame_if_appropriate();
}

extern "C" G_MODULE_EXPORT
void
on_dynamic_validation_include_missing_sidechains_checkbutton_toggled(GtkCheckButton *button,
                                                                     gpointer  user_data) {

   std::cout << "on_dynamic_validation_include_missing_sidechains_button_toggled() " << std::endl;

}



extern "C" G_MODULE_EXPORT
void
on_go_to_ligand_button_clicked(GtkButton *button,
                               gpointer   user_data) {
  go_to_ligand();
}


extern "C" G_MODULE_EXPORT
void
on_graphics_grab_focus_button_clicked (GtkButton       *button,
                                       gpointer         user_data) {
   graphics_info_t g;
   g.graphics_grab_focus();

}

#include "cc-interface-graphics.hh"

extern "C" G_MODULE_EXPORT
void
on_coot_points_button_clicked(GtkButton       *button,
                              gpointer         user_data) {
   show_coot_points_frame();
}



extern "C" G_MODULE_EXPORT
void
on_gaussian_surface_cancel_button_clicked(GtkButton       *button,
                                          gpointer         user_data) {
   GtkWidget *frame = widget_from_builder("gaussian_surface_frame");
   gtk_widget_set_visible(frame, FALSE);
}




extern "C" G_MODULE_EXPORT
void
on_gaussian_surface_ok_button_clicked(GtkButton       *button,
                                      gpointer         user_data) {

   std::cout << "read the gui - make the surface" << std::endl;
   GtkWidget *frame = widget_from_builder("gaussian_surface_frame");
   GtkWidget *mol_chooser_combobox = widget_from_builder("gaussian_surface_molecule_chooser_combobox");
   GtkWidget *e_sigma          = widget_from_builder("gaussian_surface_sigma_entry");
   GtkWidget *e_radius         = widget_from_builder("gaussian_surface_radius_entry");
   GtkWidget *e_contour_level  = widget_from_builder("gaussian_surface_contour_level_entry");
   GtkWidget *e_b_factor       = widget_from_builder("gaussian_surface_b_factor_entry");
   GtkWidget *e_chain_col_mode = widget_from_builder("gaussian_surface_chain_colour_entry");

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(mol_chooser_combobox));

   try {
      float sigma  = coot::util::string_to_float(gtk_editable_get_text(GTK_EDITABLE(e_sigma)));
      float radius = coot::util::string_to_float(gtk_editable_get_text(GTK_EDITABLE(e_radius)));
      float cl     = coot::util::string_to_float(gtk_editable_get_text(GTK_EDITABLE(e_contour_level)));
      float bf     = coot::util::string_to_float(gtk_editable_get_text(GTK_EDITABLE(e_b_factor)));
      int cc_mode  = coot::util::string_to_int(gtk_editable_get_text(GTK_EDITABLE(e_chain_col_mode)));
      set_gaussian_surface_sigma(sigma);
      set_gaussian_surface_box_radius(radius);
      set_gaussian_surface_contour_level(cl);
      set_gaussian_surface_fft_b_factor(bf);
      set_gaussian_surface_chain_colour_mode(cc_mode);
      gaussian_surface(imol);
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }

   gtk_widget_set_visible(frame, FALSE);

}

void
fill_comboboxtext_with_atom_of_residue_type(const char *rn, GtkWidget *comboboxtext) {

   std::cout << "fill_comboboxtext_with_atom_of_residue_type() --- start --- " << std::endl;

   if (rn) {
      int imol = 0;
      std::string residue_type(rn);
      auto &geom = *graphics_info_t::Geom_p();
      bool state = geom.have_dictionary_for_residue_type(residue_type, imol, graphics_info_t::cif_dictionary_read_number, true);
      graphics_info_t::cif_dictionary_read_number++;
      std::pair<bool, coot::dictionary_residue_restraints_t> rp = geom.get_monomer_restraints(residue_type, imol);
      if (rp.first) {
	 // don't do it like this -- maybe more code is needed. With this in place the combobox is just blank entries.
         // gtk_cell_layout_clear(GTK_CELL_LAYOUT(comboboxtext)); // clear combobox
         const auto &restraints = rp.second;
         const auto &atoms = restraints.atom_info;
         for (const auto &atom : atoms) {
            if (atom.type_symbol != "H") {
               std::cout << "  type " << rn << " fill comboboxtext " << comboboxtext << " with " << atom.atom_id << std::endl;
               gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(comboboxtext), atom.atom_id.c_str(), atom.atom_id.c_str());
            }
         }
      } else {
         std::cout << "ERROR get_monomer_restraints() failed for type " << residue_type << std::endl;
      }
   } else {
      std::cout << "ERROR:: null rn in fill_comboboxtext_with_atom_of_residue_type()" << std::endl;
   }
}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_first_delete_atom_checkbutton_toggled(GtkCheckButton *checkbutton,
						     gpointer         user_data) {

   std::cout << "delete atom first toggled" << std::endl;
   GtkWidget *delete_atom_combobox = widget_from_builder("acedrg_link_first_delete_atom_chooser_combobox");
   GtkWidget *entry                = widget_from_builder("acedrg_link_first_residue_name_entry");
   if (delete_atom_combobox) {
      if (entry) {
         const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (t) {
            if (gtk_check_button_get_active(checkbutton)) {
               gtk_widget_set_sensitive(delete_atom_combobox, TRUE);
               fill_comboboxtext_with_atom_of_residue_type(t, delete_atom_combobox);
            } else {
               gtk_widget_set_sensitive(delete_atom_combobox, FALSE);
            }
         }
      }
   }

}


extern "C" G_MODULE_EXPORT
void
on_acedrg_link_first_change_bond_order_checkbutton_toggled(GtkCheckButton *checkbutton,
							   gpointer         user_data) {

   std::cout << "cbo first toggled" << std::endl;
   GtkWidget *combobox = widget_from_builder("acedrg_link_first_change_bond_order_combobox");
   GtkWidget *combobox_atom_1 = widget_from_builder("acedrg_link_first_change_bond_order_atom_1_combobox");
   GtkWidget *combobox_atom_2 = widget_from_builder("acedrg_link_first_change_bond_order_atom_2_combobox");
   GtkWidget *entry           = widget_from_builder("acedrg_link_first_residue_name_entry");
   if (checkbutton) {
      if (gtk_check_button_get_active(checkbutton)) {
         gtk_widget_set_sensitive(combobox, TRUE);
         gtk_widget_set_sensitive(combobox_atom_1, TRUE);
         gtk_widget_set_sensitive(combobox_atom_2, TRUE);
         if (entry) {
            const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
            if (t) {
               fill_comboboxtext_with_atom_of_residue_type(t, combobox_atom_1);
               fill_comboboxtext_with_atom_of_residue_type(t, combobox_atom_2);
            }
         }
      } else {
         gtk_widget_set_sensitive(combobox, FALSE);
         gtk_widget_set_sensitive(combobox_atom_1, FALSE);
         gtk_widget_set_sensitive(combobox_atom_2, FALSE);
      }
   }

}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_second_delete_atom_checkbutton_toggled(GtkCheckButton *checkbutton,
						     gpointer         user_data) {

   std::cout << "delete atom second toggled" << std::endl;
   GtkWidget *delete_atom_combobox = widget_from_builder("acedrg_link_second_delete_atom_chooser_combobox");
   GtkWidget *entry                = widget_from_builder("acedrg_link_second_residue_name_entry");
   if (delete_atom_combobox) {
      if (entry) {
         const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (t) {
            if (gtk_check_button_get_active(checkbutton)) {
               gtk_widget_set_sensitive(delete_atom_combobox, TRUE);
               fill_comboboxtext_with_atom_of_residue_type(t, delete_atom_combobox);
            } else {
               gtk_widget_set_sensitive(delete_atom_combobox, FALSE);
            }
         }
      }
   }
}


extern "C" G_MODULE_EXPORT
void
on_acedrg_link_second_change_bond_order_checkbutton_toggled(GtkCheckButton *checkbutton,
							   gpointer         user_data) {

   std::cout << "cbo toggled second" << std::endl;
   GtkWidget *combobox = widget_from_builder("acedrg_link_second_change_bond_order_combobox");
   GtkWidget *combobox_atom_1 = widget_from_builder("acedrg_link_second_change_bond_order_atom_1_combobox");
   GtkWidget *combobox_atom_2 = widget_from_builder("acedrg_link_second_change_bond_order_atom_2_combobox");
   GtkWidget *entry           = widget_from_builder("acedrg_link_second_residue_name_entry");
   if (checkbutton) {
      if (gtk_check_button_get_active(checkbutton)) {
         gtk_widget_set_sensitive(combobox, TRUE);
         gtk_widget_set_sensitive(combobox_atom_1, TRUE);
         gtk_widget_set_sensitive(combobox_atom_2, TRUE);
         if (entry) {
            const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
            if (t) {
               fill_comboboxtext_with_atom_of_residue_type(t, combobox_atom_1);
               fill_comboboxtext_with_atom_of_residue_type(t, combobox_atom_2);
            }
         }
      } else {
         gtk_widget_set_sensitive(combobox, FALSE);
         gtk_widget_set_sensitive(combobox_atom_1, FALSE);
         gtk_widget_set_sensitive(combobox_atom_2, FALSE);
      }
   }
}


extern "C" G_MODULE_EXPORT
void
on_acedrg_link_first_change_charge_on_atom_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                               gpointer         user_data) {

   std::cout << "change charge atom first toggled" << std::endl;
   GtkWidget *combobox = widget_from_builder("acedrg_link_first_change_charge_on_atom_chooser_combobox");
   GtkWidget *entry    = widget_from_builder("acedrg_link_first_residue_name_entry");
   if (combobox) {
      if (entry) {
         const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (t) {
            if (gtk_check_button_get_active(checkbutton)) {
               gtk_widget_set_sensitive(combobox, TRUE);
               fill_comboboxtext_with_atom_of_residue_type(t, combobox);
            } else {
               gtk_widget_set_sensitive(combobox, FALSE);
            }
         }
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_second_change_charge_on_atom_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                                gpointer         user_data) {

   std::cout << "change charge atom second toggled" << std::endl;
   GtkWidget *combobox = widget_from_builder("acedrg_link_second_change_charge_on_atom_chooser_combobox");
   GtkWidget *entry    = widget_from_builder("acedrg_link_second_residue_name_entry");
   if (combobox) {
      if (entry) {
         const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (t) {
            if (gtk_check_button_get_active(checkbutton)) {
               gtk_widget_set_sensitive(combobox, TRUE);
               fill_comboboxtext_with_atom_of_residue_type(t, combobox);
            } else {
               gtk_widget_set_sensitive(combobox, FALSE);
            }
         }
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_ok_button_clicked(GtkButton       *button,
                                 gpointer         user_data) {

   auto simple_link = [] (const std::string &residue_name_first,
                          const std::string &residue_name_second,
                          const std::string &atom_name_first,
                          const std::string &atom_name_second,
                          const std::string &cif_file_name_1,
                          const std::string &cif_file_name_2,
                          const std::string &bond_order) {

      std::string ss = "LINK: ";

      ss += "RES-NAME-1 ";
      ss += residue_name_first;
      ss += " ";
      ss += "ATOM-NAME-1 ";
      ss += atom_name_first;
      ss += " ";
      if (!cif_file_name_1.empty())
         ss + "FILE-1 " + cif_file_name_1;

      ss += "RES-NAME-2 ";
      ss += residue_name_second;
      ss += " ";
      ss += "ATOM-NAME-2 ";
      ss += atom_name_second;
      if (!cif_file_name_2.empty())
         ss + "FILE-2 " + cif_file_name_2;

      ss += " BOND-TYPE ";
      ss += coot::util::upcase(bond_order);
      std::cout << ss << std::endl;
      run_acedrg_link_generation(ss);
   };

   // Here's an exmple:
   //
   // LINK: RES-NAME-1 LYS ATOM-NAME-1 NZ RES-NAME-2 PLP ATOM-NAME-2 C4A BOND-TYPE DOUBLE DELETE ATOM O 1 DELETE ATOM O4A 2 CHANGE BOND C OXT double 1 CHANGE BOND C4 C4A triple 2


   auto link = [] (const std::string &residue_name_first, const std::string &residue_name_second,
                   const std::string &atom_name_first, const std::string &atom_name_second,
                   const std::string &cif_file_name_1, const std::string &cif_file_name_2,
                   const std::string &bond_order,
                   bool delete_atom_first, const char *da_first,
                   bool delete_atom_second, const char *da_second,
                   bool change_bond_order_first, const char *cbo_first,
                   const char *change_bond_order_first_atom_1,  const char *change_bond_order_first_atom_2,
                   bool change_bond_order_second, const char *cbo_second,
                   const char *change_bond_order_second_atom_1, const char *change_bond_order_second_atom_2,
                   bool change_charge_on_first_residue_atom,  const char *change_charge_on_first_atom,
                   bool change_charge_on_second_residue_atom, const char *change_charge_on_second_atom) {

      std::string ss = "LINK: ";
      ss += "RES-NAME-1 ";
      ss += residue_name_first;
      ss += " ";
      ss += "ATOM-NAME-1 ";
      ss += atom_name_first;
      ss += " ";
      if (!cif_file_name_1.empty())
         ss + "FILE-1 " + cif_file_name_1 + std::string(" ");
      if (delete_atom_first)
         if (da_first)
         ss += "DELETE " + std::string(da_first) + std::string(" 1 ");
      if (change_charge_on_first_residue_atom)
         ss += std::string("CHANGE CHARGE ") + std::string(change_charge_on_first_atom) + std::string(" 1 ");
      if (change_bond_order_first)
         if (cbo_first)
            if (change_bond_order_first_atom_1)
               if (change_bond_order_first_atom_1)
                  ss += std::string("CHANGE BOND ") + std::string(change_bond_order_first_atom_1) + std::string(" ") +
                     std::string(change_bond_order_first_atom_2) + std::string(" ") + std::string(cbo_first) + " 1 ";

      ss += "RES-NAME-2 ";
      ss += residue_name_second;
      ss += " ";
      ss += "ATOM-NAME-2 ";
      ss += atom_name_second;
      if (!cif_file_name_2.empty())
         ss += std::string("FILE-2 ") + cif_file_name_2;
      if (delete_atom_second)
         if (da_second)
         ss += "DELETE " + std::string(da_second) + std::string(" 2 ");
      if (change_charge_on_second_residue_atom)
         ss += std::string("CHANGE CHARGE ") + std::string(change_charge_on_first_atom) + std::string(" 2 ");
      if (change_bond_order_second)
         if (cbo_second)
            if (change_bond_order_second_atom_1)
               if (change_bond_order_second_atom_1)
                  ss += std::string("CHANGE BOND ") + std::string(change_bond_order_second_atom_1) + std::string(" ") +
                  change_bond_order_second_atom_2 + std::string(" ") + std::string(cbo_second) + " 2 ";

      ss += std::string(" BOND-TYPE ");
      ss += coot::util::upcase(bond_order);
      run_acedrg_link_generation(ss);
      std::cout << ss << std::endl;
   };

   auto get_cif_file_name = [] (const std::string &residue_type) {

      std::string file_name;
      auto &geom = *graphics_info_t::Geom_p();
      int imol = 0;
      bool state = geom.have_dictionary_for_residue_type(residue_type, imol, graphics_info_t::cif_dictionary_read_number, true);
      graphics_info_t::cif_dictionary_read_number++;
      std::pair<bool, coot::dictionary_residue_restraints_t> rp = geom.get_monomer_restraints(residue_type, imol);
      if (rp.first) {
          file_name = rp.second.cif_file_name;
      }
      return file_name;
   };

   GtkWidget *w = widget_from_builder("acedrg_link_interface_frame");
   gtk_widget_set_visible(w, FALSE);

   GtkWidget *bond_order_combobox = widget_from_builder("acedrg_link_bond_order_combobox");

   GtkWidget *entry_first                         = widget_from_builder("acedrg_link_first_residue_name_entry");
   GtkWidget *atom_name_combobox_first            = widget_from_builder("acedrg_link_first_atom_name_chooser_combobox");

   GtkWidget *entry_second                         = widget_from_builder("acedrg_link_second_residue_name_entry");
   GtkWidget *atom_name_combobox_second            = widget_from_builder("acedrg_link_second_atom_name_chooser_combobox");

   // delete atom
   GtkWidget *delete_atom_checkbutton_second       = widget_from_builder("acedrg_link_second_delete_atom_checkbutton");
   GtkWidget *delete_atom_combobox_second          = widget_from_builder("acedrg_link_second_delete_atom_chooser_combobox");
   GtkWidget *delete_atom_checkbutton_first       = widget_from_builder("acedrg_link_first_delete_atom_checkbutton");
   GtkWidget *delete_atom_combobox_first          = widget_from_builder("acedrg_link_first_delete_atom_chooser_combobox");

   // change bond order
   GtkWidget *change_bond_order_checkbutton_first = widget_from_builder("acedrg_link_first_change_bond_order_checkbutton");
   GtkWidget *change_bond_order_combobox_first    = widget_from_builder("acedrg_link_first_change_bond_order_combobox");
   GtkWidget *change_bond_order_checkbutton_second = widget_from_builder("acedrg_link_second_change_bond_order_checkbutton");
   GtkWidget *change_bond_order_combobox_second    = widget_from_builder("acedrg_link_second_change_bond_order_combobox");
   GtkWidget *cbcbof1                              = widget_from_builder("acedrg_link_first_change_bond_order_atom_1_combobox");
   GtkWidget *cbcbof2                              = widget_from_builder("acedrg_link_first_change_bond_order_atom_2_combobox");
   GtkWidget *cbcbos1                              = widget_from_builder("acedrg_link_second_change_bond_order_atom_1_combobox");
   GtkWidget *cbcbos2                              = widget_from_builder("acedrg_link_second_change_bond_order_atom_2_combobox");

   // change charge
   GtkWidget *cc_combobox_first  = widget_from_builder("acedrg_link_first_change_charge_on_atom_chooser_combobox");
   GtkWidget *cc_combobox_second = widget_from_builder("acedrg_link_second_change_charge_on_atom_chooser_combobox");
   GtkWidget *change_charge_first_checkbutton  = widget_from_builder("on_acedrg_link_first_change_charge_on_atom_checkbutton");
   GtkWidget *change_charge_second_checkbutton = widget_from_builder("on_acedrg_link_second_change_charge_on_atom_checkbutton");

   // need to add a pair of comboboxes for change bond order atom names for both first and second.

   if (!bond_order_combobox) return;

   if (entry_first && atom_name_combobox_first && delete_atom_combobox_first && delete_atom_checkbutton_first) {
      std::cout << "Here A " << change_bond_order_combobox_first << " " << change_bond_order_checkbutton_first << std::endl;
      if (change_bond_order_combobox_first && change_bond_order_checkbutton_first) {
         std::cout << "Here B " << std::endl;
         if (entry_second && atom_name_combobox_second && delete_atom_combobox_second && delete_atom_checkbutton_second) {
            std::cout << "Here C " << std::endl;
         if (change_bond_order_combobox_second && change_bond_order_checkbutton_second) {
               std::cout << "Here D " << atom_name_combobox_first << " " << atom_name_combobox_second << std::endl;
               std::cout << "Here D first  " << GTK_IS_COMBO_BOX_TEXT(atom_name_combobox_first)  << std::endl;
               std::cout << "Here D second " << GTK_IS_COMBO_BOX_TEXT(atom_name_combobox_second) << std::endl;
               char *bond_order = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(bond_order_combobox));
               if (bond_order) {
                  std::cout << "Here D 0" << std::endl;
                  char *atom_name_first  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(atom_name_combobox_first));
                  std::cout << "Here D 1" << std::endl;
                  char *atom_name_second = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(atom_name_combobox_second));
                  std::cout << "Here D 2" << std::endl;
                  if (atom_name_first  == 0) std::cout << "Here D with atom_name_first null"  << std::endl;
                  if (atom_name_second == 0) std::cout << "Here D with atom_name_second null" << std::endl;
                  std::cout << "Here Done D" << std::endl;
                  if (atom_name_first && atom_name_second) {
                     std::cout << "Here E atom_name_first " << atom_name_first << std::endl;
                     std::cout << "Here E atom_name_second " << atom_name_second << std::endl;
                     bool change_bond_order_first  = false;
                     bool change_bond_order_second = false;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(change_bond_order_checkbutton_first)))  change_bond_order_first  = true;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(change_bond_order_checkbutton_second))) change_bond_order_second = true;
                     bool delete_atom_first  = false;
                     bool delete_atom_second = false;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(delete_atom_checkbutton_first)))  delete_atom_first  = true;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(delete_atom_checkbutton_second))) delete_atom_second = true;
                     const char *residue_name_first  = gtk_editable_get_text(GTK_EDITABLE(entry_first));
                     const char *residue_name_second = gtk_editable_get_text(GTK_EDITABLE(entry_second));
                     bool change_charge_first  = false;
                     bool change_charge_second = false;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(change_charge_first_checkbutton)))  change_charge_first  = true;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(change_charge_second_checkbutton))) change_charge_second = true;
                     if (residue_name_first) {
                        if (residue_name_second) {
                           std::string cif_file_name_1 = get_cif_file_name(residue_name_first);
                           std::string cif_file_name_2 = get_cif_file_name(residue_name_second);
                           char *da_first   = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(delete_atom_combobox_first));
                           char *da_second  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(delete_atom_combobox_second));
                           char *cbo_first  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(change_bond_order_combobox_first));
                           char *cbo_second = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(change_bond_order_combobox_second));
                           std::cout << "Here F " << change_bond_order_first << " " << change_bond_order_second << " "
                                     << delete_atom_first << " " << delete_atom_second << std::endl;
                           if (atom_name_first && atom_name_second) {
                              if (!change_bond_order_first && !change_bond_order_second && !delete_atom_first && !delete_atom_second) {
                                 simple_link(residue_name_first, residue_name_second, atom_name_first, atom_name_second,
                                             cif_file_name_1, cif_file_name_2, std::string(bond_order));
                              } else {
                                 if (cbcbof1 && cbcbof2 && cbcbos1 && cbcbos2) {
                                    const char *change_bond_order_first_atom_1  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cbcbof1));
                                    const char *change_bond_order_first_atom_2  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cbcbof2));
                                    const char *change_bond_order_second_atom_1 = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cbcbos1));
                                    const char *change_bond_order_second_atom_2 = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cbcbos2));
                                    const char *cc_first_atom                   = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cc_combobox_first));
                                    const char *cc_second_atom                  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cc_combobox_second));
                                    link(residue_name_first, residue_name_second, atom_name_first, atom_name_second,
                                         cif_file_name_1, cif_file_name_2, std::string(bond_order),
                                         delete_atom_first, da_first,
                                         delete_atom_second, da_second,
                                         change_bond_order_first, cbo_first,   change_bond_order_first_atom_1,  change_bond_order_first_atom_2,
                                         change_bond_order_second, cbo_second, change_bond_order_second_atom_1, change_bond_order_second_atom_2,
                                         change_charge_first, cc_first_atom,
                                         change_charge_second, cc_second_atom);
                                 } else {
                                    std::cout << "combobox cbo lookup failure" << std::endl;
                                 }
                              }
                           }
                        }
                     }
                  } else {
                     std::cout << "WARNING:: BAD input: atom_name_first && atom_name_second failed" << std::endl;
                  }
               }
            }
         }
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_first_residue_activate(GtkEntry *entry, gpointer user_data) {

   GtkWidget *atom_name_combobox_first = widget_from_builder("acedrg_link_first_atom_name_chooser_combobox");
   std::cout << "Fill this: " << atom_name_combobox_first << std::endl;
   const char *rn = gtk_editable_get_text(GTK_EDITABLE(entry));
   gtk_widget_set_sensitive(atom_name_combobox_first, TRUE);
   if (rn) {
      fill_comboboxtext_with_atom_of_residue_type(rn, atom_name_combobox_first);
   } else {
      std::cout << "No residue name first " << std::endl;
   }
}


extern "C" G_MODULE_EXPORT
void
on_acedrg_link_second_residue_activate(GtkEntry *entry, gpointer user_data) {

   GtkWidget *atom_name_combobox_second = widget_from_builder("acedrg_link_second_atom_name_chooser_combobox");
   std::cout << "Fill this: " << atom_name_combobox_second << std::endl;
   const char *rn = gtk_editable_get_text(GTK_EDITABLE(entry));
   gtk_widget_set_sensitive(atom_name_combobox_second, TRUE);
   if (rn) {
      fill_comboboxtext_with_atom_of_residue_type(rn, atom_name_combobox_second);
   } else {
      std::cout << "No residue name second" << std::endl;
   }
}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_cancel_button_clicked(GtkButton       *button,
				     gpointer         user_data) {

  std::cout << "Cancel" << std::endl;
  GtkWidget *w = widget_from_builder("acedrg_link_interface_frame");
  gtk_widget_set_visible(w, FALSE);
}
