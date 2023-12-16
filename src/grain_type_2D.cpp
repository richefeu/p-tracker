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

#include "grain_type_2D.hpp"

grain_type_2D::grain_type_2D() {}

// Copy Ctor
grain_type_2D::grain_type_2D(const grain_type_2D &c) {
  refcoord_xpix = c.refcoord_xpix;
  refcoord_ypix = c.refcoord_ypix;
  refrot = c.refrot;
  radius_pix = c.radius_pix;

  dx = c.dx;
  dy = c.dy;
  drot = c.drot;

  dx_prev = c.dx_prev;
  dy_prev = c.dy_prev;
  drot_prev = c.drot_prev;

  upix = c.upix;
  vpix = c.vpix;
  rot_inc = c.rot_inc;
  NCC = c.NCC;
  NCC_rescue = c.NCC_rescue;
  NCC_subpix = c.NCC_subpix;

  mean0 = c.mean0;
  C0C0 = c.C0C0;
  masked = c.masked;
  num_neighbour = c.num_neighbour;

  // other data are not copied
}

void grain_type_2D::backup() {
  dx_prev = dx;
  dy_prev = dy;
  drot_prev = drot;
}

void grain_type_2D::reset() { NCC = NCC_rescue = NCC_subpix = 0.0; }
