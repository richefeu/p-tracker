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

#include "p-tracker.hpp"

ptracker::ptracker() {
  imageData.resize(2);

  disto_parameters.resize(14);
  disto_parameters[0] = -1.0;
  disto_parameters[1] = -1.0;
  disto_parameters[2] = 0.0;
  disto_parameters[3] = 0.0;
  disto_parameters[4] = 0.0;
  disto_parameters[5] = 0.0;
  disto_parameters[6] = 0.0;
  disto_parameters[7] = 0.0;
  disto_parameters[8] = -1.0;
  disto_parameters[9] = -1.0;
  disto_parameters[10] = 0.0;
  disto_parameters[11] = 0.0;
  disto_parameters[12] = 0.0;
  disto_parameters[13] = 0.0;

  disto_parameters_perturb.resize(14);
  disto_parameters_perturb[0] = 1e-6;
  disto_parameters_perturb[1] = 1e-6;
  disto_parameters_perturb[2] = 1e-6;
  disto_parameters_perturb[3] = 1e-6;
  disto_parameters_perturb[4] = 1e-6;
  disto_parameters_perturb[5] = 1e-6;
  disto_parameters_perturb[6] = 1e-6;
  disto_parameters_perturb[7] = 1e-6;
  disto_parameters_perturb[8] = 1e-6;
  disto_parameters_perturb[9] = 1e-6;
  disto_parameters_perturb[10] = 1e-6;
  disto_parameters_perturb[11] = 1e-6;
  disto_parameters_perturb[12] = 1e-6;
  disto_parameters_perturb[13] = 1e-6;

  DemosaicModelName[0] = "LibRaw linear interpolation";
  DemosaicModelName[1] = "LibRaw VNG interpolation";
  DemosaicModelName[2] = "LibRaw PPG interpolation";
  DemosaicModelName[3] = "LibRaw AHD interpolation";
  DemosaicModelName[4] = "LibRaw DCB interpolation";
  DemosaicModelName[5] = "LibRaw demosaic pack GPL2 NOT SUPPORTED ANYMORE -> LibRaw AHD interpolation";
  DemosaicModelName[6] = "LibRaw demosaic pack GPL2 NOT SUPPORTED ANYMORE -> LibRaw AHD interpolation";
  DemosaicModelName[7] = "LibRaw demosaic pack GPL2 NOT SUPPORTED ANYMORE -> LibRaw AHD interpolation";
  DemosaicModelName[8] = "LibRaw demosaic pack GPL2 NOT SUPPORTED ANYMORE -> LibRaw AHD interpolation";
  DemosaicModelName[9] = "LibRaw demosaic pack GPL2 NOT SUPPORTED ANYMORE -> LibRaw AHD interpolation";
  DemosaicModelName[10] = "LibRaw demosaic pack GPL3 NOT SUPPORTED ANYMORE -> LibRaw AHD interpolation";
  DemosaicModelName[11] = "LibRaw DHT interpolation";
  DemosaicModelName[12] = "LibRaw Modified AHD interpolation";
  DemosaicModelName[-1] = "p-tracker averaging over 2x2 array";
  DemosaicModelName[-2] = "p-tracker weight-averaging over 3x3 array";
  DemosaicModelName[-20] = "p-tracker achromatic, no demosaising";
}

// This is the main function to read an image.
// The image can be RAW or TIFF so that the library libraw or libtiff is used.
// The image format is not recognized, but the user needs to set RawImage=1 to say that RAW format has to be used
void ptracker::read_image(int i, int num, bool first_time) {
  char name[256];
  sprintf(name, image_name, num);
  if (RawImages) {
    read_raw_image(i, name, first_time);
  } else {
    read_tiff_image(i, name, first_time);
  }
}

// Function that uses libtiff to read an image, and to compute and store gray levels in a c-style table.
void ptracker::read_tiff_image(int i, const char* name, bool first_time) {
  if (i < 0 || i > 1) {
    std::cerr << "@read_tiff_image, i = " << i << std::endl;
  }
  double tbeg = get_time();
  std::cout << "Read image named " << name << std::endl;

  // http://41j.com/blog/2011/10/simple-libtiff-example/
  TIFF* tif = TIFFOpen(name, "r");

  if (!tif) {
    std::cerr << "@read_tiff_image, cannot read tiff file named '" << name << "'" << std::endl;
    std::exit(EXIT_SUCCESS);
  }
  uint32_t w, h;
  size_t npixels;
  uint32_t* raster;

  char* infobuf;
  if (TIFFGetField(tif, TIFFTAG_DATETIME, &infobuf))
    imageData[i].dateTime = std::string(infobuf);
  else
    imageData[i].dateTime = "dateTime unknown";

  // The image data is actually not used for TIFF images
  imageData[i].iso_speed = 0.0;
  imageData[i].shutter = 0.0;
  imageData[i].aperture = 0.0;
  imageData[i].focal_len = 0.0;
  imageData[i].shot_order = 0;

  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
  std::cout << "TIFF image size: " << w << "x" << h << std::endl;
  npixels = w * h;

  raster = (uint32_t*)_TIFFmalloc(npixels * sizeof(uint32_t));
  if (raster != NULL) {
    if (!TIFFReadRGBAImage(tif, w, h, raster, 0)) {
      std::cerr << "@read_tiff_image, cannot read tiff file named '" << name << "'" << std::endl;
      std::exit(EXIT_SUCCESS);
    }
  }

  if (first_time) {
    dimx = (int)w;
    dimy = (int)h;

    // Reserve memory for image
    if (!image.empty()) {
      image.clear();
    }
    image.resize(2);
    for (size_t i = 0; i < 2; i++) {
      image[i].resize(dimx);
      for (int ix = 0; ix < dimx; ix++) {
        image[i][ix].resize(dimy);
      }
    }
  }

  double r, g, b;
  double fact = 1.0 / 255.0;
  for (int row = 0; row < dimy; ++row) {
    int shift = row * w;
    for (int col = 0; col < dimx; ++col) {
      int p = shift + col;
      r = TIFFGetR(raster[p]) * fact;
      g = TIFFGetG(raster[p]) * fact;
      b = TIFFGetB(raster[p]) * fact;
      if (r < 0.0) {
        r = 0.0;
        std::cerr << "@read_tiff_image, r < 0.0\n";
      }
      if (r > 1.0) {
        r = 1.0;
        std::cerr << "@read_tiff_image, r > 1.0\n";
      }
      if (g < 0.0) {
        g = 0.0;
        std::cerr << "@read_tiff_image, g < 0.0\n";
      }
      if (g > 1.0) {
        g = 1.0;
        std::cerr << "@read_tiff_image, g > 1.0\n";
      }
      if (b < 0.0) {
        b = 0.0;
        std::cerr << "@read_tiff_image, b < 0.0\n";
      }
      if (b > 1.0) {
        b = 1.0;
        std::cerr << "@read_tiff_image, b > 1.0\n";
      }
      // r, g and b are now values in the range [0, 1]

      // Store in a c-style table, and convert in graylevel
      // Remember that in tiff format the coordinate (0, 0) is in lower-left corner.
      image[i][col][h - row - 1] = (uint16_t)(65535.0 * (0.299 * r + 0.587 * g + 0.114 * b));
    }
  }

  _TIFFfree(raster);
  TIFFClose(tif);

  if (first_time) {
    int j = 1 - i;
    for (int y = 0; y < dimy; ++y) {
      for (int x = 0; x < dimx; ++x) {
        image[j][x][y] = image[i][x][y];
      }
    }
  }

  std::cout << "[DONE in " << get_time() - tbeg << " seconds]\n";
}

// il faut voir ici http://www.libraw.org/node/555 pour voir comment proceder
// ici aussi :
// http://stackoverflow.com/questions/22355491/libraw-is-making-my-images-too-bright-compared-to-nikons-own-converter
void ptracker::read_raw_image(int i, const char* name, bool first_time) {
  if (i < 0 || i > 1) {
    std::cerr << "@read_raw_image, i = " << i << std::endl;
  }

  LibRaw iProcessor;
  int IO_error = iProcessor.open_file(name);
  if (IO_error != 0) {
    std::cout << "\nSorry but your image named " << name << " does not exist..." << std::endl;
    std::cout << "   ..PROGRAM STOPPED.." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  double tbeg = get_time();
  std::cout << "Read and demosaicing image named " << name << std::endl;
  std::cout << "    DemosaicModel = " << DemosaicModel;
  if (DemosaicModelName.find(DemosaicModel) != DemosaicModelName.end()) {
    std::cout << " (" << DemosaicModelName[DemosaicModel] << ")";
  }

  std::cout << std::endl;

  // Get data about the image shot
  imageData[i].iso_speed = iProcessor.imgdata.other.iso_speed;
  imageData[i].shutter = iProcessor.imgdata.other.shutter;
  imageData[i].aperture = iProcessor.imgdata.other.aperture;
  imageData[i].focal_len = iProcessor.imgdata.other.focal_len;
  imageData[i].shot_order = iProcessor.imgdata.other.shot_order;
  imageData[i].dateTime = timestamp2string(iProcessor.imgdata.other.timestamp);

  std::cout << " Camera manufacturer = " << iProcessor.imgdata.idata.make << std::endl;
  std::cout << "        Camera model = " << iProcessor.imgdata.idata.model << std::endl;
  std::cout << "           iso_speed = " << imageData[i].iso_speed << std::endl;
  std::cout << "             shutter = " << imageData[i].shutter << std::endl;
  std::cout << "            aperture = " << imageData[i].aperture << std::endl;
  std::cout << "           focal_len = " << imageData[i].focal_len << std::endl;
  std::cout << "          shot_order = " << imageData[i].shot_order << std::endl;
  std::cout << "            dateTime = " << imageData[i].dateTime << std::endl;

  if (first_time) {
    dimx = iProcessor.imgdata.sizes.width;
    dimy = iProcessor.imgdata.sizes.height;
    std::cout << "RAW image size = " << dimx << "x" << dimy << std::endl;

    // Reserve memory for image
    if (!image.empty()) {
      image.clear();
    }
    image.resize(2);
    for (size_t im = 0; im < 2; im++) {
      image[im].resize(dimx);
      for (int ix = 0; ix < dimx; ix++) {
        image[im][ix].resize(dimy);
      }
    }
  }

  iProcessor.unpack();
  uint16_t MinGray = 65535, MaxGray = 0;
  double fact = 1.0 / (double)(iProcessor.imgdata.color.maximum);

  if (DemosaicModel >= 0) {  // That are those of libRaw
    iProcessor.imgdata.params.user_qual = DemosaicModel;
    iProcessor.dcraw_process();
    double r, g, b;
    int yshift = 0;
    for (int y = 0; y < dimy; ++y) {
      yshift = y * dimx;
      for (int x = 0; x < dimx; ++x) {
        // According to http://www.libraw.org/node/1995, we can ignore g2 level
        int pos = yshift + x;
        r = iProcessor.imgdata.image[pos][0] * fact;
        g = iProcessor.imgdata.image[pos][1] * fact;
        b = iProcessor.imgdata.image[pos][2] * fact;

        // Store in a c-style table, and convert in graylevel
        image[i][x][y] = (uint16_t)(iProcessor.imgdata.color.maximum * (0.299 * r + 0.587 * g + 0.114 * b));

        if (image[i][x][y] > MaxGray) MaxGray = image[i][x][y];
        if (image[i][x][y] < MinGray) MinGray = image[i][x][y];
      }
    }
  } else if (DemosaicModel == -1) {

    int width = iProcessor.imgdata.sizes.raw_width;
    int yoffset = iProcessor.imgdata.sizes.top_margin;
    int xoffset = iProcessor.imgdata.sizes.left_margin;

    double deraster;
    for (int x = 0; x < dimx; ++x) {
      for (int y = 0; y < dimy; ++y) {
        // 11
        // 11
        deraster = 0.25 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y)]);
        deraster += 0.25 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x + 1) + width * (yoffset + y)]);
        deraster +=
            0.25 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x + 1) + width * (yoffset + y + 1)]);
        deraster += 0.25 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y + 1)]);

        image[i][x][y] = (unsigned int)floor(deraster);
        if (image[i][x][y] > MaxGray) MaxGray = image[i][x][y];
        if (image[i][x][y] < MinGray) MinGray = image[i][x][y];
      }
    }
  } else if (DemosaicModel == -2) {

    int width = iProcessor.imgdata.sizes.raw_width;
    int yoffset = iProcessor.imgdata.sizes.top_margin;
    int xoffset = iProcessor.imgdata.sizes.left_margin;

    double deraster;
    for (int x = 0; x < dimx; ++x) {
      for (int y = 0; y < dimy; ++y) {
        // 121
        // 242
        // 121
        deraster =
            0.0625 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x - 1) + width * (yoffset + y - 1)]);
        deraster += 0.125 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y - 1)]);
        deraster +=
            0.0625 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x + 1) + width * (yoffset + y - 1)]);

        deraster += 0.125 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x - 1) + width * (yoffset + y)]);
        deraster += 0.25 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y)]);
        deraster += 0.125 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x + 1) + width * (yoffset + y)]);

        deraster +=
            0.0625 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x - 1) + width * (yoffset + y + 1)]);
        deraster += 0.125 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y + 1)]);
        deraster +=
            0.0625 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x + 1) + width * (yoffset + y + 1)]);

        image[i][x][y] = (unsigned int)floor(deraster);
        if (image[i][x][y] > MaxGray) MaxGray = image[i][x][y];
        if (image[i][x][y] < MinGray) MinGray = image[i][x][y];
      }
    }
  } else if (DemosaicModel == -20) {  // for achromatic raw images (no demosaicing)

    int width = iProcessor.imgdata.sizes.raw_width;
    int yoffset = iProcessor.imgdata.sizes.top_margin;
    int xoffset = iProcessor.imgdata.sizes.left_margin;

    for (int x = 0; x < dimx; ++x) {
      for (int y = 0; y < dimy; ++y) {
        image[i][x][y] = iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y)];
        if (image[i][x][y] > MaxGray) MaxGray = image[i][x][y];
        if (image[i][x][y] < MinGray) MinGray = image[i][x][y];
      }
    }
  }

  iProcessor.recycle();

  std::cout << "MinGray value = " << MinGray << '\n';
  std::cout << "MaxGray value = " << MaxGray << '\n';

  // rescale according to MaxGray and MinGray.
  // It should be ok since we use ZNCC (which is not sensitive to a scale factor)
  // The default behaviour is to NOT rescale the gray levels
  if (rescaleGrayLevels) {
    std::cout << "Image greylevels are rescaled !!!\n";
    double fact = 65535.0 / (double)(MaxGray - MinGray);
    for (int y = 0; y < dimy; ++y) {
      for (int x = 0; x < dimx; ++x) {
        image[i][x][y] = (uint16_t)floor((image[i][x][y] - MinGray) * fact);
      }
    }
  }

  if (first_time) {
    int j = 1 - i;
    for (int y = 0; y < dimy; ++y) {
      for (int x = 0; x < dimx; ++x) {
        image[j][x][y] = image[i][x][y];
      }
    }
  }

  std::cout << "[DONE in " << get_time() - tbeg << " seconds, user_qual = " << iProcessor.imgdata.params.user_qual
            << "]\n";
}

void ptracker::header() {
  // clang-format off
  std::cout << '\n';
  std::cout << "    _/_/_/_/       _/_/_/_/_/  _/_/_/      _/_/      _/_/_/  _/    _/  _/_/_/_/  _/_/_/ \n";
  std::cout << "   _/      _/         _/      _/    _/  _/    _/  _/        _/  _/    _/        _/    _/\n";
  std::cout << "  _/_/_/_/    _/_/   _/      _/_/_/    _/_/_/_/  _/        _/_/      _/_/_/    _/_/_/   \n";
  std::cout << " _/                 _/      _/    _/  _/    _/  _/        _/  _/    _/        _/    _/  \n";
  std::cout << "_/                 _/      _/    _/  _/    _/    _/_/_/  _/    _/  _/_/_/_/  _/    _/   \n";
  std::cout << std::endl;
  // clang-format on
  std::cout << "   p-tracker  Copyright (C) 2023  Vincent Richefeu, Gael Combe\n";
  std::cout << "   This program comes with ABSOLUTELY NO WARRANTY\n";
  // std::cout << " vincent.richefeu@univ-grenoble-alpes.fr, gael.combe@3sr-grenoble.fr " << std::endl;
  std::cout << std::endl;
  std::cout << "   Version number: " << P_TRACKER_VERSION << std::endl;

#if defined(_OPENMP)
  std::cout << "   Multi-threading: enabled (OpenMP)" << std::endl;
#else
  std::cout << "   Multi-threading: Disabled" << std::endl;
#endif

  std::cout << std::endl;
}

double ptracker::get_time() {
#if defined(_OPENMP)
  return omp_get_wtime();
#else
  return (double)std::clock() / (double)CLOCKS_PER_SEC;
#endif
}

// The particle tracking procedure (default procedure)
void ptracker::particle_tracking() {
  std::cout << "*** PROCEDURE PARTICLE IMAGE TRACKING ***" << std::endl;

  read_image(im_index_ref, iref, true);  // Read reference image

  double tbeg;
  int igrain;
  int irescue;

  if (use_neighbour_list) {
    std::cout << "Build neighbour list ... " << std::flush;
    tbeg = get_time();
    for (int igrain = 0; igrain < num_grains; igrain++) {
      grain[igrain].neighbour.clear();
      find_neighbours(igrain);
    }
    std::cout << "[DONE in " << get_time() - tbeg << " seconds]" << std::endl;
  }

  for (int num_image = ibeg; num_image <= iend; num_image += iinc) {

    std::cout << "\n____ Correlations from image " << iref << " towards image " << num_image << " ("
              << num_image - ibeg + 1 << "/" << iend - ibeg + 1 << ")\n";
    read_image(im_index_current, num_image);

    // Re-build the neighbour list periodically
    if (use_neighbour_list && num_image % period_rebuild_neighbour_list == 0) {
      std::cout << "Rebuilding the neighbour list ... " << std::flush;
      tbeg = get_time();
      for (int igrain = 0; igrain < num_grains; igrain++) {
        find_neighbours(igrain);
      }

      std::cout << "[DONE in " << get_time() - tbeg << " seconds]" << std::endl;
    }

    // Precomputations
    if (require_precomputations) {
      std::cout << "Precomputations ... " << std::flush;
      tbeg = get_time();
      do_precomputations();
      std::cout << "[DONE in " << get_time() - tbeg << " seconds]" << std::endl;
    }

    // coarse pixel correlation (with rotations).
    // The very first attempts
    num_to_be_rescued = num_to_be_super_rescued = 0;
    std::cout << "Follow " << num_grains << " grains ... " << std::endl;

    if (rescue_level >= 0) {
      tbeg = get_time();

      progress = 0;
#pragma omp parallel for schedule(dynamic)
      for (igrain = 0; igrain < num_grains; igrain++) {
        grain[igrain].reset();   // reset les NCC
        grain[igrain].backup();  // sauvegarde les dx, dy et drot
        if (grain[igrain].masked) continue;
        follow_pattern_pixel(igrain);
#pragma omp critical
        { loadbar(++progress, grain.size()); }
      }

      std::cerr << std::endl;

      std::cout << "[DONE in " << get_time() - tbeg << " seconds]" << std::endl;
    }

    // According to a minimum value of NCC, a first rescue is attempted
    if (rescue_level >= 1) {
      if (num_to_be_rescued > 0) std::cout << "Try to rescue " << num_to_be_rescued << " grains" << std::endl;

      tbeg = get_time();

      progress = 0;
#pragma omp parallel for private(igrain) schedule(dynamic)
      for (irescue = 0; irescue < num_to_be_rescued; irescue++) {
        igrain = to_be_rescued[irescue];
        follow_pattern_rescue_pixel(igrain);
#pragma omp critical
        { loadbar(++progress, num_to_be_rescued); }
      }
      std::cerr << std::endl;

      // Display
      for (irescue = 0; irescue < num_to_be_rescued; irescue++) {
        igrain = to_be_rescued[irescue];
        fprintf(stdout, "\nGrain %6d \t NCC = %6.4f  ->  ", igrain, grain[igrain].NCC);
        fprintf(stdout, "%6.4f \t diff = %6.4f \t ", grain[igrain].NCC_rescue,
                grain[igrain].NCC_rescue - grain[igrain].NCC);
        if (grain[igrain].NCC_rescue > NCC_min)
          fprintf(stdout, "[\033[32mOK\033[0m]\n");
        else
          fprintf(stdout, "[\033[31mFAIL\033[0m]\n");
      }

      if (num_to_be_rescued > 0) std::cout << "[rescue DONE in " << get_time() - tbeg << " seconds]" << std::endl;
    }

    // According to a minimum value of NCC, a last super_rescue is attempted
    if (rescue_level >= 2) {
      if (num_to_be_super_rescued > 0)
        std::cout << "Try to super-rescue " << num_to_be_super_rescued << " grains" << std::endl;

      tbeg = get_time();

      progress = 0;
#pragma omp parallel for private(igrain) schedule(dynamic)
      for (irescue = 0; irescue < num_to_be_super_rescued; irescue++) {
        igrain = to_be_super_rescued[irescue];
        follow_pattern_super_rescue_pixel(igrain);
#pragma omp critical
        { loadbar(++progress, num_to_be_rescued); }
      }
      std::cerr << std::endl;

      // Display
      for (irescue = 0; irescue < num_to_be_super_rescued; irescue++) {
        igrain = to_be_super_rescued[irescue];
        fprintf(stdout, "Grain %6d \t NCC = %6.4f  ->  ", igrain, grain[igrain].NCC);
        fprintf(stdout, "%6.4f \t diff = %6.4f \t ", grain[igrain].NCC_rescue,
                grain[igrain].NCC_rescue - grain[igrain].NCC);
        if (grain[igrain].NCC_rescue > NCC_min_super)
          fprintf(stdout, "[\033[32mOK\033[0m]\n");
        else {
          fprintf(stdout, "[\033[31mFAIL\033[0m]\n");
          fprintf(stdout, " Please, fix the pb manually\n");
        }
      }

      if (num_to_be_super_rescued > 0)
        std::cout << "[super-rescue DONE in " << get_time() - tbeg << " seconds]" << std::endl;
    }

    // Sub-pixel
    if (subpixel) {
      std::cout << "Sub-pixel resolution ..." << std::endl;

      tbeg = get_time();

      progress = 0;
#pragma omp parallel for schedule(dynamic)
      for (igrain = 0; igrain < num_grains; igrain++) {
        if (grain[igrain].masked) {
          continue;
        }
        follow_pattern_subpixel_xyR(igrain);
#pragma omp critical
        { loadbar(++progress, num_grains); }
      }
      std::cerr << std::endl;
      std::cout << "[Sub-pixel resolution DONE in " << get_time() - tbeg << " seconds]" << std::endl;
    }

    for (igrain = 0; igrain < num_grains; igrain++) {
      grain[igrain].upix = grain[igrain].dx - grain[igrain].dx_prev;
      grain[igrain].vpix = grain[igrain].dy - grain[igrain].dy_prev;
      grain[igrain].rot_inc = grain[igrain].drot - grain[igrain].drot_prev;
    }

    // Save the DIC file
    save_grains(num_image);

    // Change reference (Not yet tested!!!)
    if (num_image - iref >= iraz) {  // equal in fact!
      std::cout << '\n';
      std::cout << "***********************************************\n\n";
      std::cout << " Changing reference image from " << iref << " to " << num_image << "\n\n";
      std::cout << "***********************************************" << std::endl;

      int tmp = im_index_current;
      im_index_current = im_index_ref;
      im_index_ref = tmp;  // Swap the first index for the array 'image'
      double actual_pos;
      for (igrain = 0; igrain < num_grains; igrain++) {
        actual_pos = (double)(grain[igrain].refcoord_xpix) + grain[igrain].dx;
        grain[igrain].refcoord_xpix = std::round(actual_pos);
        grain[igrain].dx = actual_pos - grain[igrain].refcoord_xpix;

        actual_pos = (double)(grain[igrain].refcoord_ypix) + grain[igrain].dy;
        grain[igrain].refcoord_ypix = std::round(actual_pos);
        grain[igrain].dy = actual_pos - grain[igrain].refcoord_ypix;

        grain[igrain].refrot += grain[igrain].drot;
        grain[igrain].drot = 0.0;
      }
      iref = num_image;
      require_precomputations = true;  // so that precomputations are updated
    }
  }
}

void ptracker::rotate_pixel_pattern(int igrain, int i, double c, double s, int* xpixel, int* ypixel) {
  // Note that the y is oriented downwards, which is why the trigonometric rotation is the opposite
  double xpix = c * grain[igrain].pattern[i].dx + s * grain[igrain].pattern[i].dy;
  double ypix = -s * grain[igrain].pattern[i].dx + c * grain[igrain].pattern[i].dy;
  *xpixel = std::round(xpix);
  *ypixel = std::round(ypix);
}

// Pre-computation of mean0, C0C0 and interpolated_zone for each grain
void ptracker::do_precomputations() {
  double c0, s0;
  int xpixel0, ypixel0;
  double diffC0;

  for (int igrain = 0; igrain < num_grains; igrain++) {
    int refcoordx = grain[igrain].refcoord_xpix;
    int refcoordy = grain[igrain].refcoord_ypix;
    double refrot = grain[igrain].refrot;

    grain[igrain].mean0 = 0.0;
    grain[igrain].C0C0 = 0.0;
    c0 = cos(refrot);
    s0 = sin(refrot);
    for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
      rotate_pixel_pattern(igrain, i, c0, s0, &xpixel0, &ypixel0);
      grain[igrain].pattern0_rotated[i].dx = xpixel0;
      grain[igrain].pattern0_rotated[i].dy = ypixel0;
      xpixel0 += refcoordx;
      ypixel0 += refcoordy;
      grain[igrain].mean0 += (double)image[im_index_ref][xpixel0][ypixel0];
    }
    grain[igrain].mean0 /= (double)grain[igrain].pattern.size();
    for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
      xpixel0 = grain[igrain].pattern0_rotated[i].dx + refcoordx;
      ypixel0 = grain[igrain].pattern0_rotated[i].dy + refcoordy;
      diffC0 = (double)image[im_index_ref][xpixel0][ypixel0] - grain[igrain].mean0;
      grain[igrain].C0C0 += diffC0 * diffC0;
    }
  }  // Loop over igrain

  require_precomputations = false;
}

// Based on Highest NCC
void ptracker::follow_pattern_pixel(int igrain) {
  int refcoordx = grain[igrain].refcoord_xpix;
  int refcoordy = grain[igrain].refcoord_ypix;
  double refrot = grain[igrain].refrot;
  double c1, s1;
  int xpixel0, ypixel0;

  int grain_dx = std::round(grain[igrain].dx);
  int grain_dy = std::round(grain[igrain].dy);
  double rest_dx = grain[igrain].dx - grain_dx;
  double rest_dy = grain[igrain].dy - grain_dy;
  grain[igrain].upix = grain[igrain].vpix = grain[igrain].rot_inc = 0.0;

  double drot;     // Increment test pour l'angle de rotation
  int upix, vpix;  // Increment test pour la translation

  double mean1, C0C1, C1C1;

  int xpixel1, ypixel1;
  double total_rot;

  double NCC_test, best_NCC = 0.0;
  double best_drot = 0.0;
  int best_upix = 0.0, best_vpix = 0.0;
  double diffC1;
  double rotmax = fabs(search_zone.inc_rot * search_zone.num_rot);

  for (drot = -rotmax; drot <= rotmax; drot += search_zone.inc_rot) {
    total_rot = refrot + grain[igrain].drot + drot;
    c1 = cos(total_rot);
    s1 = sin(total_rot);
    for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
      rotate_pixel_pattern(igrain, i, c1, s1, &xpixel1, &ypixel1);
      grain[igrain].pattern1_rotated[i].dx = xpixel1;
      grain[igrain].pattern1_rotated[i].dy = ypixel1;
    }
    for (upix = -search_zone.left; upix <= search_zone.right; upix++) {
      for (vpix = -search_zone.up; vpix <= search_zone.down; vpix++) {

        mean1 = 0.0;
        C0C1 = C1C1 = 0.0;
        for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
          xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
          ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
          mean1 += (double)image[im_index_current][xpixel1][ypixel1];
        }
        mean1 /= (double)grain[igrain].pattern.size();
        for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
          xpixel0 = grain[igrain].pattern0_rotated[i].dx + refcoordx;
          ypixel0 = grain[igrain].pattern0_rotated[i].dy + refcoordy;
          xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
          ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
          diffC1 = (double)image[im_index_current][xpixel1][ypixel1] - mean1;
          C0C1 += ((double)image[im_index_ref][xpixel0][ypixel0] - grain[igrain].mean0) * diffC1;
          C1C1 += diffC1 * diffC1;
        }

        NCC_test = C0C1 / sqrt(grain[igrain].C0C0 * C1C1);

        if (best_NCC < NCC_test) {
          best_NCC = NCC_test;
          best_upix = (double)upix;
          best_vpix = (double)vpix;
          best_drot = drot;
        }

      }  // Loop vpix
    }    // Loop upix
  }      // Loop rot

  grain[igrain].NCC = best_NCC;
  grain[igrain].NCC_rescue = best_NCC;

  if (best_NCC < NCC_min) {
#pragma omp critical
    { to_be_rescued[num_to_be_rescued++] = igrain; }
  } else {
    grain[igrain].upix = best_upix - rest_dx;
    grain[igrain].vpix = best_vpix - rest_dy;
    grain[igrain].rot_inc = best_drot;
    grain[igrain].dx = grain_dx + best_upix;
    grain[igrain].dy = grain_dy + best_vpix;
    grain[igrain].drot += best_drot;
  }
}

void ptracker::follow_pattern_rescue_pixel(int igrain) {
  int refcoordx = grain[igrain].refcoord_xpix;
  int refcoordy = grain[igrain].refcoord_ypix;
  double refrot = grain[igrain].refrot;
  double c, s;
  int xpixel0, ypixel0;

  int grain_dx = std::round(grain[igrain].dx);
  int grain_dy = std::round(grain[igrain].dy);
  double rest_dx = grain[igrain].dx - grain_dx;
  double rest_dy = grain[igrain].dy - grain_dy;

  double drot;     // Increment test pour l'angle de rotation
  int upix, vpix;  // Increment test pour la translation

  double mean1, C0C1, C1C1;

  int i_allowed;
  bool is_allowed;

  int xpixel1, ypixel1;
  double total_rot;

  double NCC_test, best_NCC = 0.0;
  double best_drot = 0.0;
  int best_upix = 0.0, best_vpix = 0.0;
  double rotmax = fabs(search_zone_rescue.inc_rot * search_zone_rescue.num_rot);

  // We define a list of relative coordinates used to search igrain in the image 1.
  // To do that, the list of neighbour of igrain is used
  const int nbmax =
      (search_zone_rescue.right + search_zone_rescue.left + 1) * (search_zone_rescue.up + search_zone_rescue.down + 1);
  std::vector<relative_coord_type> rescue_allowed(nbmax);

  int radi_i = (int)grain[igrain].radius_pix;

  int nb_allowed = 0;

  if (use_neighbour_list) {
    for (int ix = -search_zone_rescue.left; ix <= search_zone_rescue.right; ix++) {
      for (int iy = -search_zone_rescue.up; iy <= search_zone_rescue.down; iy++) {
        is_allowed = true;
        int xcentre_i = refcoordx + grain_dx + ix;
        int ycentre_i = refcoordy + grain_dy + iy;
        for (int i_neigh = 0; i_neigh < grain[igrain].num_neighbour; i_neigh++) {
          int jgrain = grain[igrain].neighbour[i_neigh];
          double jgrain_NCC = grain[jgrain].NCC_rescue;
          if (jgrain_NCC > NCC_min) {
            int xcentre_j = grain[jgrain].refcoord_xpix + std::round(grain[jgrain].dx);
            int ycentre_j = grain[jgrain].refcoord_ypix + std::round(grain[jgrain].dy);
            int radi_j = (int)grain[jgrain].radius_pix;
            double dstx = xcentre_j - xcentre_i;
            double dsty = ycentre_j - ycentre_i;
            double dist2 = dstx * dstx + dsty * dsty;
            double sum_radii = (radi_i + radi_j) * 0.9;  // only 90% of the radius declared
            double sum_radii2 = sum_radii * sum_radii;
            if (dist2 < sum_radii2) {
              is_allowed = false;
              break;
            }
          }
        }
        if (is_allowed) {
          rescue_allowed[nb_allowed].dx = ix;
          rescue_allowed[nb_allowed].dy = iy;
          nb_allowed++;
        }
      }
    }
  } else {
    for (int ix = -search_zone_rescue.left; ix <= search_zone_rescue.right; ix++) {
      for (int iy = -search_zone_rescue.up; iy <= search_zone_rescue.down; iy++) {
        rescue_allowed[nb_allowed].dx = ix;
        rescue_allowed[nb_allowed].dy = iy;
        nb_allowed++;
      }
    }
  }

  for (drot = -rotmax; drot <= rotmax; drot += search_zone_rescue.inc_rot) {
    total_rot = refrot + grain[igrain].drot + drot;
    c = cos(total_rot);
    s = sin(total_rot);
    for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
      rotate_pixel_pattern(igrain, i, c, s, &xpixel1, &ypixel1);
      grain[igrain].pattern1_rotated[i].dx = xpixel1;
      grain[igrain].pattern1_rotated[i].dy = ypixel1;
    }

    for (i_allowed = 0; i_allowed < nb_allowed; i_allowed++) {
      upix = rescue_allowed[i_allowed].dx;
      vpix = rescue_allowed[i_allowed].dy;
      mean1 = 0.0;
      C0C1 = C1C1 = 0.0;
      for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
        xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
        ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
        mean1 += image[im_index_current][xpixel1][ypixel1];
      }
      mean1 /= (double)grain[igrain].pattern.size();
      for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
        xpixel0 = grain[igrain].pattern0_rotated[i].dx + refcoordx;
        ypixel0 = grain[igrain].pattern0_rotated[i].dy + refcoordy;
        xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
        ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
        C0C1 += (image[im_index_ref][xpixel0][ypixel0] - grain[igrain].mean0) *
                (image[im_index_current][xpixel1][ypixel1] - mean1);
        C1C1 +=
            (image[im_index_current][xpixel1][ypixel1] - mean1) * (image[im_index_current][xpixel1][ypixel1] - mean1);
      }

      NCC_test = C0C1 / sqrt(grain[igrain].C0C0 * C1C1);

      if (best_NCC < NCC_test) {
        best_NCC = NCC_test;
        best_upix = (double)upix;
        best_vpix = (double)vpix;
        best_drot = drot;
      }

    }  // Loop i_allowed
  }    // Loop rot

  grain[igrain].NCC_rescue = best_NCC;

  if (best_NCC < NCC_min_super) {
#pragma omp critical
    { to_be_super_rescued[num_to_be_super_rescued++] = igrain; }
  } else {
    grain[igrain].upix = best_upix - rest_dx;
    grain[igrain].vpix = best_vpix - rest_dy;
    grain[igrain].rot_inc = best_drot;
    grain[igrain].dx = grain_dx + best_upix;
    grain[igrain].dy = grain_dy + best_vpix;
    grain[igrain].drot += best_drot;
  }
}

void ptracker::follow_pattern_super_rescue_pixel(int igrain) {
  int refcoordx = grain[igrain].refcoord_xpix;
  int refcoordy = grain[igrain].refcoord_ypix;
  double refrot = grain[igrain].refrot;
  double c, s;
  int xpixel0, ypixel0;

  int grain_dx = std::round(grain[igrain].dx);
  int grain_dy = std::round(grain[igrain].dy);
  double rest_dx = grain[igrain].dx - grain_dx;
  double rest_dy = grain[igrain].dy - grain_dy;

  double drot;     // Increment test pour l'angle de rotation
  int upix, vpix;  // Increment test pour la translation
  double mean1, C0C1, C1C1;

  int i_allowed;
  bool is_allowed;

  int xpixel1, ypixel1;
  double total_rot;

  double NCC_test, best_NCC = 0.0;
  double best_drot = 0.0;
  int best_upix = 0.0, best_vpix = 0.0;
  double rotmax = fabs(search_zone_super_rescue.inc_rot * search_zone_super_rescue.num_rot);

  // We define a list of relative coordinates used to search igrain in the image 1.
  // To do that, the list of neighbour of igrain is used
  const int nbmax = (search_zone_super_rescue.right + search_zone_super_rescue.left + 1) *
                    (search_zone_super_rescue.up + search_zone_super_rescue.down + 1);
  std::vector<relative_coord_type> super_rescue_allowed(nbmax);

  int radi_i = (int)grain[igrain].radius_pix;

  int nb_allowed = 0;

  if (use_neighbour_list) {
    for (int ix = -search_zone_super_rescue.left; ix <= search_zone_super_rescue.right; ix++) {
      for (int iy = -search_zone_super_rescue.up; iy <= search_zone_super_rescue.down; iy++) {
        is_allowed = true;
        int xcentre_i = refcoordx + grain_dx + ix;
        int ycentre_i = refcoordy + grain_dy + iy;
        for (int i_neigh = 0; i_neigh < grain[igrain].num_neighbour; i_neigh++) {
          int jgrain = grain[igrain].neighbour[i_neigh];
          double jgrain_NCC = std::max(grain[jgrain].NCC, grain[jgrain].NCC_rescue);
          if (jgrain_NCC > NCC_min_super) {
            int xcentre_j = grain[jgrain].refcoord_xpix + std::round(grain[jgrain].dx);
            int ycentre_j = grain[jgrain].refcoord_ypix + std::round(grain[jgrain].dy);
            int radi_j = (int)grain[jgrain].radius_pix;
            double dstx = xcentre_j - xcentre_i;
            double dsty = ycentre_j - ycentre_i;
            double dist2 = dstx * dstx + dsty * dsty;
            double sum_radii = (radi_i + radi_j) * 0.9;  // only 90% of the radius declared
            double sum_radii2 = sum_radii * sum_radii;
            if (dist2 < sum_radii2) {
              is_allowed = false;
              break;
            }
          }
        }
        if (is_allowed) {
          super_rescue_allowed[nb_allowed].dx = ix;
          super_rescue_allowed[nb_allowed].dy = iy;
          nb_allowed++;
        }
      }
    }
  } else {
    for (int ix = -search_zone_super_rescue.left; ix <= search_zone_super_rescue.right; ix++) {
      for (int iy = -search_zone_super_rescue.up; iy <= search_zone_super_rescue.down; iy++) {
        super_rescue_allowed[nb_allowed].dx = ix;
        super_rescue_allowed[nb_allowed].dy = iy;
        nb_allowed++;
      }
    }
  }

  for (drot = -rotmax; drot <= rotmax; drot += search_zone_super_rescue.inc_rot) {
    total_rot = refrot + grain[igrain].drot + drot;
    c = cos(total_rot);
    s = sin(total_rot);
    for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
      rotate_pixel_pattern(igrain, i, c, s, &xpixel1, &ypixel1);
      grain[igrain].pattern1_rotated[i].dx = xpixel1;
      grain[igrain].pattern1_rotated[i].dy = ypixel1;
    }

    for (i_allowed = 0; i_allowed < nb_allowed; i_allowed++) {
      upix = super_rescue_allowed[i_allowed].dx;
      vpix = super_rescue_allowed[i_allowed].dy;
      mean1 = 0.0;
      C0C1 = C1C1 = 0.0;
      for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
        xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
        ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
        mean1 += (double)image[im_index_current][xpixel1][ypixel1];
      }
      mean1 /= (double)grain[igrain].pattern.size();
      for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
        xpixel0 = grain[igrain].pattern0_rotated[i].dx + refcoordx;
        ypixel0 = grain[igrain].pattern0_rotated[i].dy + refcoordy;
        xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
        ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
        C0C1 += ((double)image[im_index_ref][xpixel0][ypixel0] - grain[igrain].mean0) *
                ((double)image[im_index_current][xpixel1][ypixel1] - mean1);
        C1C1 += ((double)image[im_index_current][xpixel1][ypixel1] - mean1) *
                ((double)image[im_index_current][xpixel1][ypixel1] - mean1);
      }

      NCC_test = C0C1 / sqrt(grain[igrain].C0C0 * C1C1);

      if (best_NCC < NCC_test) {
        best_NCC = NCC_test;
        best_upix = (double)upix;
        best_vpix = (double)vpix;
        best_drot = drot;
      }
    }

  }  // Loop for drot

  grain[igrain].NCC_rescue = best_NCC;

  if (best_NCC > NCC_min_super) {  // If the super rescue has no effect, then we don't change the movement increment
    grain[igrain].upix = best_upix - rest_dx;
    grain[igrain].vpix = best_vpix - rest_dy;
    grain[igrain].rot_inc = best_drot;

    grain[igrain].dx = grain_dx + best_upix;
    grain[igrain].dy = grain_dy + best_vpix;
    grain[igrain].drot += best_drot;
  }
}

double ptracker::NCC_to_minimize_xyR(std::vector<double>& X) {

#if defined(_OPENMP)
  int thread_id = omp_get_thread_num();
#else
  int thread_id = 0;
#endif

  int igrain = igrain_of_thread[thread_id];

  double inv_num_pix = 1.0 / (double)(grain[igrain].pattern.size());
  int xpix0 = grain[igrain].refcoord_xpix;
  int ypix0 = grain[igrain].refcoord_ypix;
  double mean0 = grain[igrain].mean0;
  double C0C0 = grain[igrain].C0C0;

  double xpix1 = (double)grain[igrain].refcoord_xpix + (double)grain[igrain].dx + X[0];
  double ypix1 = (double)grain[igrain].refcoord_ypix + (double)grain[igrain].dy + X[1];

  double rot = grain[igrain].drot + X[2];  // from refrot which can be non-zero
  double c1 = cos(rot), s1 = sin(rot);

  std::vector<double> interpol_values;

  // Compute mean values
  double xsubpixel1, ysubpixel1;
  double mean1 = 0.0, int_grey_val;
  double x, y;

  for (size_t i = 0; i < grain[igrain].pattern0_rotated.size(); i++) {
    x = grain[igrain].pattern0_rotated[i].dx;
    y = grain[igrain].pattern0_rotated[i].dy;

    // Note that if y is pointing downwards, the trigonometric sens of rotation must be reversed.
    xsubpixel1 = xpix1 + c1 * (double)x + s1 * (double)y;
    ysubpixel1 = ypix1 - s1 * (double)x + c1 * (double)y;

    int_grey_val = IMAGE_INTERPOLATOR->getValue(image[im_index_current], xsubpixel1, ysubpixel1);
    interpol_values.push_back(int_grey_val);
    mean1 += int_grey_val;
  }
  mean1 *= inv_num_pix;

  // NCC computation
  double C0C1 = 0.0, C1C1 = 0.0;
  double diffC0, diffC1;
  int xpixel0, ypixel0;
  int num_interpol_values = 0;
  for (size_t i = 0; i < grain[igrain].pattern0_rotated.size(); i++) {
    x = grain[igrain].pattern0_rotated[i].dx;
    y = grain[igrain].pattern0_rotated[i].dy;

    diffC1 = interpol_values[num_interpol_values++] - mean1;
    xpixel0 = xpix0 + x;
    ypixel0 = ypix0 + y;

    diffC0 = (double)image[im_index_ref][xpixel0][ypixel0] - mean0;
    C0C1 += diffC0 * diffC1;
    C1C1 += diffC1 * diffC1;
  }

  double NCC = C0C1 / sqrt(C0C0 * C1C1);

  double W = 1.0;
  double subpix_displacement = sqrt(X[0] * X[0] + X[1] * X[1]);
  if (subpix_displacement > subpix_displacement_max) {
    W = overflow_stiffness * (subpix_displacement - subpix_displacement_max) + 1.0;
  }

  return (1.0 - NCC) * W;
}

void ptracker::follow_pattern_subpixel_xyR(int igrain) {
  std::vector<double> X(3);   // The parameters to be optimized
  std::vector<double> DX(3);  // The initial variation of the params in the optimisation process

  // === Initial guess of the subpixel part of displacement and rotation
  X[0] = X[1] = X[2] = 0.0;

  // === Initial direction (to search minimum in powell method)
  DX[0] = initial_direction_dx;
  DX[1] = initial_direction_dy;
  DX[2] = initial_direction_dr;

  // === Preparation of the minimization
#if defined(_OPENMP)
  int thread_id = omp_get_thread_num();
#else
  int thread_id = 0;
#endif

  igrain_of_thread[thread_id] = igrain;

  // === Minimization
  NCC_optimizer func;
  func.ptrackerObject = this;
  Powell<NCC_optimizer> powell(func, subpix_tol);

  X = powell.minimize(X, DX);

  // === Update
  grain[igrain].NCC_subpix = 1.0 - powell.fret;

  grain[igrain].upix += X[0];
  grain[igrain].vpix += X[1];
  grain[igrain].rot_inc += X[2];

  grain[igrain].dx += X[0];
  grain[igrain].dy += X[1];
  grain[igrain].drot += X[2];
}

void ptracker::find_neighbours(int igrain) {
  double squared_dist_pix = neighbour_dist_pix * neighbour_dist_pix;

  double x_igrain = grain[igrain].refcoord_xpix + grain[igrain].dx;
  double y_igrain = grain[igrain].refcoord_ypix + grain[igrain].dy;

  grain[igrain].num_neighbour = 0;
  for (int jgrain = 0; jgrain < num_grains; jgrain++) {
    if (jgrain == igrain) {
      continue;
    }
    double x_jgrain = grain[jgrain].refcoord_xpix + grain[jgrain].dx;
    double y_jgrain = grain[jgrain].refcoord_ypix + grain[jgrain].dy;
    double dist2 = (x_jgrain - x_igrain) * (x_jgrain - x_igrain) + (y_jgrain - y_igrain) * (y_jgrain - y_igrain);
    if (dist2 <= squared_dist_pix) {
      grain[igrain].neighbour.push_back(jgrain);
    }
  }
}

int ptracker::init(int argc, char* argv[]) {

  if (argc > 2) {
    fprintf(stderr, "Usage: %s command_file\n", argv[0]);
    fprintf(stderr, "Type %s -h for help\n", argv[0]);
    std::exit(EXIT_SUCCESS);
  }

  if (argc == 2) {
    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
      header();
      std::cout << std::endl;
      std::cout << "Help can be found in the folder 'examples' or 'doc'" << std::endl;
      std::cout << "  -v, --version   Display information about the soft" << std::endl;
      std::cout << "  -h, --help      This help" << std::endl;
      std::cout << "  -s, --sizes     Print size of common types (in bits) for this computer" << std::endl;
      std::cout << std::endl;
      std::exit(EXIT_SUCCESS);
    } else if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0) {
      header();
      std::exit(EXIT_SUCCESS);
    } else if (strcmp(argv[1], "-s") == 0 || strcmp(argv[1], "--sizes") == 0) {
      std::cout << std::endl;
      std::cout << "unsigned char..... " << sizeof(unsigned char) * 8 << " bits" << std::endl;
      std::cout << "unsigned short.... " << sizeof(unsigned short) * 8 << " bits" << std::endl;
      std::cout << "unsigned int...... " << sizeof(unsigned int) * 8 << " bits" << std::endl;
      std::cout << "float............. " << sizeof(float) * 8 << " bits" << std::endl;
      std::cout << "double............ " << sizeof(double) * 8 << " bits" << std::endl;
      std::cout << "uint16_t.......... " << sizeof(uint16_t) * 8 << " bits" << std::endl;
      std::cout << "uint16_t.......... " << sizeof(uint16_t) * 8 << " bits" << std::endl;
      std::cout << std::endl;
      std::exit(EXIT_SUCCESS);
    } else {
      header();
      read_data(argv[1]);
    }
  } else {  // tracker has been invocked by double-clics or without arguments
    header();
    read_data("commands.txt");
  }

  if (rotations == 0) {  // to be replaced by if (search_zone_super_rescue.inc_rot == 0.0 || ...) (FIXME)
    std::cout << "Rotation are NOT TRACKED!" << std::endl;
    // Il faut choisir une valeur non nulle pour inc_rot sinon -> boucle infinie
    search_zone_super_rescue.inc_rot = 1.0;
    search_zone_rescue.inc_rot = 1.0;
    search_zone.inc_rot = 1.0;

    search_zone_super_rescue.num_rot = 0;
    search_zone_rescue.num_rot = 0;
    search_zone.num_rot = 0;
  }

#ifdef _OPENMP
  omp_set_num_threads(wanted_num_threads);
#endif

  std::cout << "Number of threads: " << wanted_num_threads << std::endl;

  return 1;
}

// Lecture des positions des grains "tracked"
int ptracker::read_grains(const char* name, bool ptracker_format) {
  std::ifstream grain_file(name);
  if (!grain_file) {
    std::cerr << "@ptracker::read_grains, cannot open file " << name << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  grain_file >> num_grains;

  if (!grain.empty()) {
    grain.clear();
  }
  grain.resize(num_grains);
  if (!to_be_rescued.empty()) {
    to_be_rescued.clear();
  }
  to_be_rescued.resize(num_grains);
  if (!to_be_super_rescued.empty()) {
    to_be_super_rescued.clear();
  }
  to_be_super_rescued.resize(num_grains);

  if (ptracker_format == true) {  // File with the same format as 'dic_out_x.txt'
    for (int i = 0; i < num_grains; i++) {
      // clang-format off
      grain_file >> grain[i].refcoord_xpix  // 1
                 >> grain[i].refcoord_ypix  // 2
                 >> grain[i].refrot         // 3
                 >> grain[i].radius_pix     // 4
                 >> grain[i].dx             // 5
                 >> grain[i].dy             // 6
                 >> grain[i].drot           // 7
                 >> grain[i].upix           // 8
                 >> grain[i].vpix           // 9
                 >> grain[i].rot_inc        // 10
                 >> grain[i].NCC            // 11
                 >> grain[i].NCC_rescue     // 12
                 >> grain[i].NCC_subpix;    // 13
      // clang-format on
    }
  } else {  // Format plus simple
    double radius_pix_real;
    double xpix_real, ypix_real;
    for (int i = 0; i < num_grains; i++) {
      // clang-format off
      grain_file >> xpix_real        // 1
                 >> ypix_real        // 2
                 >> grain[i].refrot  // 3
                 >> radius_pix_real; // 4
      // clang-format on
      grain[i].radius_pix = radius_pix_real;

      // Note that if the coordinates are real, then grain.dx is non-zero.
      // This adds an artificial initial displacement.
      grain[i].refcoord_xpix = (int)floor(xpix_real);
      grain[i].dx = grain[i].upix = xpix_real - floor(xpix_real);
      grain[i].refcoord_ypix = (int)floor(ypix_real);
      grain[i].dy = grain[i].vpix = ypix_real - floor(ypix_real);
    }
  }

  // Try to find the 1G2E related datasets
  // frame.try_to_plug(this);

  return 1;
}

void ptracker::make_grid(double xmin, double xmax, double ymin, double ymax, int nx, int ny, int aleaMax) {
  int dx = (int)floor((xmax - xmin) / (double)(nx - 1));
  int dy = (int)floor((ymax - ymin) / (double)(ny - 1));

  // Memory allocation
  num_grains = nx * ny;
  if (!grain.empty()) grain.clear();
  grain.resize(num_grains);
  if (!to_be_rescued.empty()) {
    to_be_rescued.clear();
  }
  to_be_rescued.resize(num_grains);
  if (!to_be_super_rescued.empty()) {
    to_be_super_rescued.clear();
  }
  to_be_super_rescued.resize(num_grains);

  std::cout << "Create grid with " << num_grains << " points" << std::endl;

  int i = 0;
  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      double angle = (double)rand() / (double)RAND_MAX * (2.0 * M_PI);  // angle in range [0 2PI]
      double alea = (double)rand() / (double)RAND_MAX * aleaMax;
      grain[i].refcoord_xpix = xmin + ix * dx + (int)std::round(alea * cos(angle));
      grain[i].refcoord_ypix = ymin + iy * dy + (int)std::round(alea * sin(angle));
      grain[i].refrot = 0.0;
      grain[i].radius_pix = 2.0;
      i++;
    }
  }
}

int ptracker::save_grains(const char* name, int num, bool simpleVersion) {
  std::ofstream grain_file_out(name);
  if (!grain_file_out) {
    std::cerr << "Cannot open file " << name << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  if (simpleVersion) {
    grain_file_out << grain.size() << std::endl;
    for (size_t i = 0; i < grain.size(); i++) {
      grain_file_out << grain[i].refcoord_xpix + grain[i].dx << ' ' << grain[i].refcoord_ypix + grain[i].dy << ' '
                     << grain[i].refrot << ' ' << grain[i].radius_pix << '\n';
    }
  } else {
    grain_file_out << num_grains << std::endl;
    for (int i = 0; i < num_grains; i++) {
      // clang-format off
      grain_file_out << std::scientific << std::setprecision(15) 
                     << std::setw(6)  << std::left  << grain[i].refcoord_xpix << ' '   // 1
                     << std::setw(6)  << std::left  << grain[i].refcoord_ypix << ' '   // 2
                     << std::setw(23) << std::right << grain[i].refrot        << ' '   // 3
                     << std::setw(23) << std::right << grain[i].radius_pix    << ' '   // 4
                     << std::setw(23) << std::right << grain[i].dx            << ' '   // 5
                     << std::setw(23) << std::right << grain[i].dy            << ' '   // 6
                     << std::setw(23) << std::right << grain[i].drot          << ' '   // 7
                     << std::setw(23) << std::right << grain[i].upix          << ' '   // 8
                     << std::setw(23) << std::right << grain[i].vpix          << ' '   // 9
                     << std::setw(23) << std::right << grain[i].rot_inc       << ' '   // 10
                     << std::setw(23) << std::right << grain[i].NCC           << ' '   // 11
                     << std::setw(23) << std::right << grain[i].NCC_rescue    << ' '   // 12
                     << std::setw(23) << std::right << grain[i].NCC_subpix    << '\n'; // 13
      // clang-format on
    }

    grain_file_out << std::endl;

    char fname[256];
    sprintf(fname, image_name, iref);
    grain_file_out << "#_FROM_image: " << fname << '\n';
    grain_file_out << "#___dateTime: " << imageData[im_index_ref].dateTime << '\n';
    sprintf(fname, image_name, num);
    grain_file_out << "#___TO_image: " << fname << '\n';
    grain_file_out << "#___dateTime: " << imageData[im_index_current].dateTime << '\n';
    grain_file_out << std::flush;
  }

  return 1;
}

int ptracker::save_grains(int num) {
  char name[256];
  sprintf(name, "dic_out_%d.txt", num);

  return save_grains(name, num);
}

void ptracker::make_rect_pattern(int igrain, int half_width_pix, int half_height_pix) {
  // Memory allocation
  int dim = (2 * half_width_pix + 1) * (2 * half_height_pix + 1);
  grain[igrain].pattern.resize(dim);
  grain[igrain].pattern0_rotated.resize(dim);
  grain[igrain].pattern1_rotated.resize(dim);

  int i = 0;
  for (int dy = -half_height_pix; dy <= half_height_pix; ++dy) {
    for (int dx = -half_width_pix; dx <= half_width_pix; ++dx) {
      grain[igrain].pattern[i].dx = dx;
      grain[igrain].pattern[i].dy = dy;
      i++;
    }
  }
}

void ptracker::make_circ_pattern(int igrain, int radius_pix) {
  // Memory allocation
  int dim = 4 * radius_pix * radius_pix;  // surdim !
  grain[igrain].pattern.resize(dim);
  grain[igrain].pattern0_rotated.resize(dim);
  grain[igrain].pattern1_rotated.resize(dim);
  int i = 0;
  double dst;
  for (int dy = -radius_pix; dy <= radius_pix; ++dy) {
    for (int dx = -radius_pix; dx <= radius_pix; ++dx) {
      dst = sqrt((double)dx * (double)dx + (double)dy * (double)dy);
      if (dst <= (double)radius_pix) {
        grain[igrain].pattern[i].dx = dx;
        grain[igrain].pattern[i].dy = dy;
        i++;
      }
    }
  }

  grain[igrain].pattern.resize(i);
  grain[igrain].pattern0_rotated.resize(i);
  grain[igrain].pattern1_rotated.resize(i);
}

void ptracker::make_ring_pattern(int igrain, int radius_IN_pix, int radius_OUT_pix) {
  // Memory allocation
  int dim = 4 * radius_OUT_pix * radius_OUT_pix;  // more than necessary!
  grain[igrain].pattern.resize(dim);
  grain[igrain].pattern0_rotated.resize(dim);
  grain[igrain].pattern1_rotated.resize(dim);

  int i = 0;
  double dst;
  for (int dy = -radius_OUT_pix; dy <= radius_OUT_pix; ++dy) {
    for (int dx = -radius_OUT_pix; dx <= radius_OUT_pix; ++dx) {
      dst = sqrt((double)dx * (double)dx + (double)dy * (double)dy);
      if (dst <= (double)radius_OUT_pix && dst >= (double)radius_IN_pix) {
        grain[igrain].pattern[i].dx = dx;
        grain[igrain].pattern[i].dy = dy;
        i++;
      }
    }
  }

  grain[igrain].pattern.resize(i);
  grain[igrain].pattern0_rotated.resize(i);
  grain[igrain].pattern1_rotated.resize(i);
}

// Format:
// number_of_patterns
// for each pattern {
// 	number_of_relative_coords
// 	for each coord {
//		xrel yrel
// 	}
// }
void ptracker::make_custom_pattern(const char* name) {
  std::ifstream pattern_file(name);
  if (!pattern_file) {
    fprintf(stderr, "Cannot open file %s\n", name);
    std::exit(EXIT_SUCCESS);
  }

  int num_patterns;
  pattern_file >> num_patterns;
  if (num_patterns != num_grains) {
    fprintf(stderr, "The number of patterns is not equal to the number of grains!\n");
    fprintf(stderr, "num_patterns : %d\n", num_patterns);
    fprintf(stderr, "num_grains : %d\n", num_grains);
    std::exit(EXIT_SUCCESS);
  }
  fprintf(stdout, "Number of patterns %d\n", num_patterns);

  int num_coords = 0;
  int dx, dy;
  for (int igrain = 0; igrain < num_patterns; igrain++) {
    pattern_file >> num_coords;
    // Reserver memoire pour pattern, pattern0_rotation et pattern1_rotation
    grain[igrain].pattern.resize(num_coords);
    grain[igrain].pattern0_rotated.resize(num_coords);
    grain[igrain].pattern1_rotated.resize(num_coords);

    // grain[igrain].num_point_pattern = num_coords;
    for (int i = 0; i < num_coords; i++) {
      pattern_file >> dx >> dy;
      grain[igrain].pattern[i].dx = dx;
      grain[igrain].pattern[i].dy = dy;
    }
  }
}

void ptracker::mask_rect(int xmin, int xmax, int ymin, int ymax) {
  for (int igrain = 0; igrain < num_grains; igrain++) {
    if (grain[igrain].refcoord_xpix >= xmin && grain[igrain].refcoord_xpix <= xmax &&
        grain[igrain].refcoord_ypix >= ymin && grain[igrain].refcoord_ypix <= ymax)
      grain[igrain].masked = true;
  }
}

// This function reads the command file
int ptracker::read_data(const char* name) {
  std::ifstream command_file(name);
  if (!command_file) {
    std::cerr << "Cannot open file " << name << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::string token;
  command_file >> token;

  while (command_file) {

    if (token[0] == '!' || token[0] == '#')
      getline(command_file, token);
    else if (token == "image_name") {
      command_file >> image_name;
    } else if (token == "dic_name") {
      command_file >> dic_name;
    } else if (token == "RawImages") {
      command_file >> RawImages;
    } else if (token == "DemosaicModel") {
      command_file >> DemosaicModel;
    } else if (token == "rescaleGrayLevels") {
      command_file >> rescaleGrayLevels;
    } else if (token == "grains_to_follow" || token == "start_file") {
      char filename[256];
      command_file >> filename;
      read_grains(filename, false);
    } else if (token == "restart_file") {
      char filename[256];
      command_file >> filename;
      read_grains(filename, true);
    } else if (token == "make_grid") {
      double xmin, xmax, ymin, ymax;
      int nx, ny;
      int aleaMax;
      command_file >> xmin >> xmax >> ymin >> ymax >> nx >> ny >> aleaMax;
      make_grid(xmin, xmax, ymin, ymax, nx, ny, aleaMax);
    } else if (token == "iref") {
      command_file >> iref;
    } else if (token == "ibeg") {
      command_file >> ibeg;
    } else if (token == "iend") {
      command_file >> iend;
    } else if (token == "iinc") {
      command_file >> iinc;
    } else if (token == "iraz") {
      command_file >> iraz;
    } else if (token == "idelta") {
      command_file >> idelta;
    }

    else if (token == "wanted_num_threads") {
      int n;
      command_file >> n;
      wanted_num_threads = fabs(n);
    } else if (token == "procedure") {
      command_file >> procedure;
    } else if (token == "verbose_level") {
      command_file >> verbose_level;
    } else if (token == "rescue_level") {
      command_file >> rescue_level;
    } else if (token == "subpixel") {
      command_file >> subpixel;
    } else if (token == "rotations") {
      command_file >> rotations;
    } else if (token == "optim_algo") {
      std::string keyword;
      command_file >> keyword;
      if (keyword == "POWELL") {
        optim_algo = POWELL;
      } else if (keyword == "POWELL_ITER") {
        optim_algo = POWELL_ITER;
      } else if (keyword == "LEVMAR") {
        optim_algo = LEVMAR;
      } else {
        std::cout << "WARNING, keyword " << keyword << " is unknown for optim_algo\n";
      }
    } else if (token == "powell_num_iter") {
      command_file >> powell_num_iter;
    }

    else if (token == "num_neighbour_max") {
      command_file >> num_neighbour_max;
    } else if (token == "neighbour_dist_pix") {
      command_file >> neighbour_dist_pix;
    } else if (token == "period_rebuild_neighbour_list") {
      command_file >> period_rebuild_neighbour_list;
    }

    else if (token == "image_interpolator") {
      command_file >> token;
      if (token == "linear") {
        IMAGE_INTERPOLATOR = &IMAGE_INTERPOLATOR_LINEAR;
      } else if (token == "cubic") {
        IMAGE_INTERPOLATOR = &IMAGE_INTERPOLATOR_CUBIC;
      } else if (token == "quintic") {
        IMAGE_INTERPOLATOR = &IMAGE_INTERPOLATOR_QUINTIC;
      }
    } else if (token == "subpix_tol") {
      command_file >> subpix_tol;
    } else if (token == "initial_direction_dx") {
      command_file >> initial_direction_dx;
    } else if (token == "initial_direction_dy") {
      command_file >> initial_direction_dy;
    } else if (token == "initial_direction_dr") {
      command_file >> initial_direction_dr;
    } else if (token == "overflow_stiffness") {
      command_file >> overflow_stiffness;
    } else if (token == "subpix_displacement_max") {
      command_file >> subpix_displacement_max;
    }

    else if (token == "search_zone.left") {
      command_file >> search_zone.left;
    } else if (token == "search_zone.right") {
      command_file >> search_zone.right;
    } else if (token == "search_zone.up") {
      command_file >> search_zone.up;
    } else if (token == "search_zone.down") {
      command_file >> search_zone.down;
    } else if (token == "search_zone.inc_rot") {
      command_file >> search_zone.inc_rot;
    } else if (token == "search_zone.num_rot") {
      command_file >> search_zone.num_rot;
    }

    else if (token == "NCC_min") {
      command_file >> NCC_min;
    } else if (token == "search_zone_rescue.left") {
      command_file >> search_zone_rescue.left;
    } else if (token == "search_zone_rescue.right") {
      command_file >> search_zone_rescue.right;
    } else if (token == "search_zone_rescue.up") {
      command_file >> search_zone_rescue.up;
    } else if (token == "search_zone_rescue.down") {
      command_file >> search_zone_rescue.down;
    } else if (token == "search_zone_rescue.inc_rot") {
      command_file >> search_zone_rescue.inc_rot;
    } else if (token == "search_zone_rescue.num_rot") {
      command_file >> search_zone_rescue.num_rot;
    }

    else if (token == "NCC_min_super") {
      command_file >> NCC_min_super;
    } else if (token == "search_zone_super_rescue.left") {
      command_file >> search_zone_super_rescue.left;
    } else if (token == "search_zone_super_rescue.right") {
      command_file >> search_zone_super_rescue.right;
    } else if (token == "search_zone_super_rescue.up") {
      command_file >> search_zone_super_rescue.up;
    } else if (token == "search_zone_super_rescue.down") {
      command_file >> search_zone_super_rescue.down;
    } else if (token == "search_zone_super_rescue.inc_rot") {
      command_file >> search_zone_super_rescue.inc_rot;
    } else if (token == "search_zone_super_rescue.num_rot") {
      command_file >> search_zone_super_rescue.num_rot;
    } else if (token == "use_neighbour_list") {
      command_file >> use_neighbour_list;
    }

    else if (token == "targetRadiusPattern") {
      command_file >> targetRadiusPattern;
    } else if (token == "pattern") {
      command_file >> token;
      if (token == "circ") {
        int radius;
        command_file >> radius;
        // Same pattern for each grain
        for (int i = 0; i < num_grains; i++)
          if (targetRadiusPattern < 0.0 || grain[i].radius_pix == targetRadiusPattern) make_circ_pattern(i, radius);
      } else if (token == "grain_diameter") {
        double radius_reduction, min_radius;
        command_file >> radius_reduction >> min_radius;
        // circular pattern for each grain with a radius dependant to grain radius
        for (int i = 0; i < num_grains; i++) {
          int radius = (int)(grain[i].radius_pix * radius_reduction);
          if (radius < min_radius) {
            radius = min_radius;
          }
          if (targetRadiusPattern < 0.0 || grain[i].radius_pix == targetRadiusPattern) make_circ_pattern(i, radius);
        }
      } else if (token == "ring") {
        int radius_IN, radius_OUT;
        command_file >> radius_IN >> radius_OUT;
        // Same pattern for each grain
        for (int i = 0; i < num_grains; i++)
          if (targetRadiusPattern < 0.0 || grain[i].radius_pix == targetRadiusPattern)
            make_ring_pattern(i, radius_IN, radius_OUT);
      } else if (token == "file") {
        char filename[256];
        command_file >> filename;
        make_custom_pattern(filename);
      } else if (token == "rect") {
        int half_width_pix, half_height_pix;
        command_file >> half_width_pix >> half_height_pix;
        // Same pattern for each grain
        for (int i = 0; i < num_grains; i++)
          if (targetRadiusPattern < 0.0 || grain[i].radius_pix == targetRadiusPattern)
            make_rect_pattern(i, half_width_pix, half_height_pix);
      }
    }

    else if (token == "mask_rect") {
      int xmin, xmax, ymin, ymax;
      command_file >> xmin >> xmax >> ymin >> ymax;
      mask_rect(xmin, xmax, ymin, ymax);
    }

    else if (token == "image_numbers_corrDisto") {
      size_t nb;
      command_file >> nb;
      image_numbers_corrDisto.resize(nb);
      for (size_t i = 0; i < nb; i++) {
        command_file >> image_numbers_corrDisto[i];
      }
    } else if (token == "imposedDisplApprox_corrDisto") {
      imposed_displ_x_pix.resize(image_numbers_corrDisto.size());
      imposed_displ_y_pix.resize(image_numbers_corrDisto.size());
      imposed_displ_x_pix[0] = imposed_displ_y_pix[0] = 0.0;
      for (size_t i = 1; i < image_numbers_corrDisto.size(); i++) {
        command_file >> imposed_displ_x_pix[i] >> imposed_displ_y_pix[i];
      }
    }

    else if (token == "equiProj_dist_range") {
      command_file >> equiProj_dist_min >> equiProj_dist_max;
    }

    else if (token == "xc_corrDistor") {
      command_file >> disto_parameters[0];
    } else if (token == "yc_corrDistor") {
      command_file >> disto_parameters[1];
    } else if (token == "K1_corrDistor") {
      command_file >> disto_parameters[2];
    } else if (token == "K2_corrDistor") {
      command_file >> disto_parameters[3];
    } else if (token == "K3_corrDistor") {
      command_file >> disto_parameters[4];
    } else if (token == "P1_corrDistor") {
      command_file >> disto_parameters[5];
    } else if (token == "P2_corrDistor") {
      command_file >> disto_parameters[6];
    } else if (token == "P3_corrDistor") {
      command_file >> disto_parameters[7];
    } else if (token == "xc_corrParallax") {
      command_file >> disto_parameters[8];
    } else if (token == "yc_corrParallax") {
      command_file >> disto_parameters[9];
    } else if (token == "angle1_corrParallax") {
      command_file >> disto_parameters[10];
    } else if (token == "A1_corrParallax") {
      command_file >> disto_parameters[11];
    } else if (token == "angle2_corrParallax") {
      command_file >> disto_parameters[12];
    } else if (token == "A2_corrParallax") {
      command_file >> disto_parameters[13];
    }

    else if (token == "dxc_corrDistor") {
      command_file >> disto_parameters_perturb[0];
    } else if (token == "dyc_corrDistor") {
      command_file >> disto_parameters_perturb[1];
    } else if (token == "dK1_corrDistor") {
      command_file >> disto_parameters_perturb[2];
    } else if (token == "dK2_corrDistor") {
      command_file >> disto_parameters_perturb[3];
    } else if (token == "dK3_corrDistor") {
      command_file >> disto_parameters_perturb[4];
    } else if (token == "dP1_corrDistor") {
      command_file >> disto_parameters_perturb[5];
    } else if (token == "dP2_corrDistor") {
      command_file >> disto_parameters_perturb[6];
    } else if (token == "dP3_corrDistor") {
      command_file >> disto_parameters_perturb[7];
    } else if (token == "dxc_corrParallax") {
      command_file >> disto_parameters_perturb[8];
    } else if (token == "dyc_corrParallax") {
      command_file >> disto_parameters_perturb[9];
    } else if (token == "dangle1_corrParallax") {
      command_file >> disto_parameters_perturb[10];
    } else if (token == "dA1_corrParallax") {
      command_file >> disto_parameters_perturb[11];
    } else if (token == "dangle2_corrParallax") {
      command_file >> disto_parameters_perturb[12];
    } else if (token == "dA2_corrParallax") {
      command_file >> disto_parameters_perturb[13];
    }

    else {
      fprintf(stdout, "Unknown token: %s\n", token.c_str());
    }

    command_file >> token;
  }

  return 1;
}

// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
void ptracker::loadbar(size_t x, size_t n, size_t w) {
  if ((x != n) && (x % (n / 100 + 1) != 0)) return;

  float ratio = x / (float)n;
  size_t c = ratio * w;

  std::cerr << std::setw(3) << (size_t)(ratio * 100) << "% |";
  for (size_t x = 0; x < c; x++) std::cerr << "|";
  for (size_t x = c; x < w; x++) std::cerr << " ";
  std::cerr << "|\r" << std::flush;
}

std::string ptracker::timestamp2string(time_t rawtime) {
  struct tm* timeinfo;
  timeinfo = localtime(&rawtime);
  char retChar[256];
  sprintf(retChar, "%s", asctime(timeinfo));
  return std::string(retChar);
}

void ptracker::run() {
  if (procedure == "particle_tracking") {
    particle_tracking();
  }
  if (procedure == "correction_distortion") {
    correction_distortion();
  } else {  // The default procedure
    particle_tracking();
  }
}
