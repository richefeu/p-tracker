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

// Research zone (for pixel-precision tracking)
struct search_zone_type {
  int left;       // Towards decreasing x (positive value)
  int right;      // Towards increasing x (positive value)
  int up;         // Towards decreasing y (positive value)
  int down;       // Towards increasing y (positive value)
  int num_rot;    // Number of tested angles (above and after, total = 2*num_rot+1)
  double inc_rot; // Increment of rotation
};


