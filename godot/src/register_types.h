/*
 * godot/src/register_types.h
 */
#ifndef COOT_GODOT_REGISTER_TYPES_H
#define COOT_GODOT_REGISTER_TYPES_H

#include <godot_cpp/core/class_db.hpp>

void initialize_coot_module(godot::ModuleInitializationLevel p_level);
void uninitialize_coot_module(godot::ModuleInitializationLevel p_level);

#endif // COOT_GODOT_REGISTER_TYPES_H
