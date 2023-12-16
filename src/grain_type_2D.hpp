// p-tracker
// Copyright (C) 2023  Vincent Richefeu, Gael Combe
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <vector>

#include "relative_coord_type.hpp"

struct grain_type_2D {
  // initial variables (generally the input data)
  int refcoord_xpix{0};   // reference x-coordinate of the point being tracked
  int refcoord_ypix{0};   // reference y-coordinate of the point being tracked
  double refrot{0.0};     // initial (reference) angle of the pattern
  double radius_pix{0.0}; // grain radius expressed in sub-pixels (input data)

  // Cumulative displacements and rotations relative to the initial variables (relative to refxxx)
  double dx{0.0};
  double dy{0.0};
  double drot{0.0};

  // Backup of cumulative movements and rotations
  double dx_prev{0.0};
  double dy_prev{0.0};
  double drot_prev{0.0};

  // Kinematic data obtain on current DIC (num_image to num_image + iinc)
  double upix{0.0};    // x translation resulting from last DIC
  double vpix{0.0};    // y translation resulting from last DIC
  double rot_inc{0.0}; // Increment of rotation resulting from last DIC

  // Correlation coefficients
  double NCC{0.0};        // Best NCC
  double NCC_rescue{0.0}; // Best NCC after rescue (0 if no rescue has been performed)
  double NCC_subpix{0.0}; // Best NCC after subpix (0 if no subpix refinement has been performed)
  
  // Precomputed data for the computation of NCCs
  double mean0{0.0};
  double C0C0{0.0};   
  
  bool masked{false}; // Pour les point dans une zone à ne pas correlé (utilisation avec une grille)

  std::vector<relative_coord_type> pattern;          // A list of relative positions that represent a pattern
  std::vector<relative_coord_type> pattern0_rotated; // Rotated pattern for the reference position
  std::vector<relative_coord_type> pattern1_rotated; // Rotated pattern for the tested position

  std::vector<int> neighbour; // List of neighbouring grains
  int num_neighbour{0};       // Number of neighbour for igrain (according to neighbour_dist_pix) @fixme supprimer

  //Ctor
  grain_type_2D();
  
  // Copy constructor
  grain_type_2D(const grain_type_2D &c);

  void backup();
  void reset();
};
