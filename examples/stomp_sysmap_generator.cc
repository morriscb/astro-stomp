#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <stomp.h>
#include <gflags/gflags.h>

namespace Stomp {
  class SysPixel : public Pixel {
  public:
    SysPixel() {
      SetResolution(0);

      SetPixnumFromXY(0, 0);
      SetWeight(0.0);

      sum_seeing_ = 0.0;
      sum_extinction_ = 0.0;
      sum_sky_ = 0.0;
      n_obj_ = 0;
    };
    SysPixel(uint32_t resolution, uint32_t pixnum, double unmasked_fraction) {
      SetResolution(resolution);

      uint32_t tmp_y = pixnum/(Stomp::Nx0*Resolution());
      uint32_t tmp_x = pixnum - Stomp::Nx0*Resolution()*tmp_y;

      SetPixnumFromXY(tmp_x, tmp_y);
      SetWeight(unmasked_fraction);

      sum_seeing_ = 0.0;
      sum_extinction_ = 0.0;
      sum_sky_ = 0.0;
      n_obj_ = 0;
    };
    SysPixel(uint32_t resolution, uint32_t pixnum, double unmasked_fraction,
	     double mean_seeing, double mean_extinction, double mean_sky,
	     int n_objects) {
      SetResolution(resolution);

      uint32_t tmp_y = pixnum/(Stomp::Nx0*Resolution());
      uint32_t tmp_x = pixnum - Stomp::Nx0*Resolution()*tmp_y;

      SetPixnumFromXY(tmp_x, tmp_y);
      SetWeight(unmasked_fraction);

      sum_seeing_ = mean_seeing*n_objects;
      sum_extinction_ = mean_extinction*n_objects;
      sum_sky_ = mean_sky*n_objects;
      n_obj_ = n_objects;
    };
    ~SysPixel() {
      SetResolution(0);

      SetPixnumFromXY(0, 0);
      SetWeight(0.0);

      sum_seeing_ = 0.0;
      sum_extinction_ = 0.0;
      sum_sky_ = 0.0;
      n_obj_ = 0;
    }
    inline void AddSysParameters(double seeing, double extinction, double sky) {
      sum_seeing_ += seeing;
      sum_extinction_ += extinction;
      sum_sky_ += sky;
      n_obj_++;
    };
    inline double MeanSeeing() {
      return (n_obj_ > 0 ? sum_seeing_/n_obj_ : -1.0);
    };
    inline double MeanExtinction() {
      return (n_obj_ > 0 ? sum_extinction_/n_obj_ : -1.0);
    };
    inline double MeanSky() {
      return (n_obj_ > 0 ? sum_sky_/n_obj_ : -1.0);
    };
    inline double UnmaskedFraction() {
      return Weight();
    };
    inline int NObjects() {
      return n_obj_;
    };

  private:
    double sum_seeing_, sum_extinction_, sum_sky_;
    int n_obj_;
  };
} // end namespace Stomp

typedef std::map<const uint32_t, Stomp::SysPixel> SysDict;
typedef SysDict::iterator SysDictIterator;

// Define our command-line flags.
DEFINE_string(input_files, "",
              "CSV list of input systematics files");
DEFINE_string(map_file, "",
              "Name of the ASCII file containing the input map.");
DEFINE_string(sysmap_file, "",
              "Name of the ASCII file to write systematics map to.");
DEFINE_int32(resolution, 256,
	     "Resolution to use for creating the systematics map.");
DEFINE_bool(single_index, false, "Use older single-index file input format.");
DEFINE_bool(no_weight, false, "Input map file is missing weight column.");

int main(int argc, char **argv) {
  std::string usage = "Usage: ";
  usage += argv[0];
  usage += " --input_files=<CSV list of input systematics files>";
  usage += " --map_file=<StompMap ASCII>";
  usage += " --sysmap_file=<SysMap ASCII>";
  google::SetUsageMessage(usage);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_input_files.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_map_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  if (FLAGS_sysmap_file.empty()) {
    std::cout << usage << "\n";
    std::cout << "Type '" << argv[0] << " --help' for a list of options.\n";
    exit(1);
  }

  // Before starting, the internal way that we're representing the systematics
  // map (with the SysDict STL map) relies on using the Stomp::Pixel Pixnum()
  // method as a key.  This method only works for resolution vales of 2048 or
  // lower, so we need to catch the user before they try using an illegal
  // value.
  if (FLAGS_resolution > 2048) {
    std::cout << "Maximum resolution is 2048.  Please lower expectations.\n";
    exit(1);
  }

  // First, we read our STOMP map into a map object.  There are a couple
  // permutations based on the various map formats that are out there: with
  // or without a weight column or in the single index or double index format.
  std::cout << "Reading in initial map...\n";
  Stomp::Map* stomp_map;
  if (FLAGS_single_index) {
    if (FLAGS_no_weight) {
      stomp_map = new Stomp::Map(FLAGS_map_file, false, false);
    } else {
      stomp_map = new Stomp::Map(FLAGS_map_file, false);
    }
  } else {
    if (FLAGS_no_weight) {
      stomp_map = new Stomp::Map(FLAGS_map_file, true, false);
    } else {
      stomp_map = new Stomp::Map(FLAGS_map_file);
    }
  }
  std::cout << "Read map from " << FLAGS_map_file << "; initial area: " <<
    stomp_map->Area() << " sq. deg.\n";

  // To prime our map, we find the coverage map of the input geometry at the
  // systematics map resolution.
  Stomp::PixelVector coverage_map;
  stomp_map->Coverage(coverage_map, FLAGS_resolution);

  stomp_map->Clear();
  delete stomp_map;

  // Now we create a map object from pixel index to SysPixel...
  SysDict sysmap;

  //... and initialize it with our coverage pixels.
  for (Stomp::PixelIterator iter=coverage_map.begin();
       iter!=coverage_map.end();++iter) {
    Stomp::SysPixel pix(iter->Resolution(), iter->Pixnum(), iter->Weight());
    sysmap[iter->Pixnum()] = pix;
  }

  // Now we start reading in our systematics files.
  std::vector<std::string> input_files;

  Stomp::Tokenize(FLAGS_input_files, input_files, ",");

  std::cout << "Parsing " << input_files.size() << " files...\n";

  uint32_t n_obj = 0, n_keep = 0;
  for (std::vector<std::string>::iterator file_iter=input_files.begin();
       file_iter!=input_files.end();++file_iter) {

    std::cout << "\tParsing " << file_iter->c_str() << "...\n";
    std::ifstream systematics_file(file_iter->c_str());
    double lambda, eta, seeing, extinction, sky;
    SysDictIterator sys_iter;

    while (!systematics_file.eof()) {
      systematics_file >> lambda >> eta >> seeing >> extinction >> sky;

      if (!systematics_file.eof()) {
	n_obj++;
	Stomp::AngularCoordinate ang(lambda, eta,
				     Stomp::AngularCoordinate::Survey);
	Stomp::Pixel pix(ang, FLAGS_resolution, 1.0);
	sys_iter = sysmap.find(pix.Pixnum());
	if (sys_iter != sysmap.end()) {
	  sys_iter->second.AddSysParameters(seeing, extinction, sky);
	  n_keep++;
	}
      }
    }
    systematics_file.close();
  }
  std::cout << "Kept " << n_keep << "/" << n_obj << " systematics objects...\n";

  std::cout << "Writing systematics map to " << FLAGS_sysmap_file << "\n";

  std::ofstream output_file(FLAGS_sysmap_file.c_str());
  // Now we iterate over our coverage pixels and write out the results
  for (Stomp::PixelIterator iter=coverage_map.begin();
       iter!=coverage_map.end();++iter) {
    SysDictIterator sys_iter = sysmap.find(iter->Pixnum());
    if (sys_iter != sysmap.end()) {
      output_file << sys_iter->second.HPixnum() << " " <<
	sys_iter->second.Superpixnum() << " " <<
	sys_iter->second.Resolution() << " " <<
	sys_iter->second.UnmaskedFraction() << " " <<
	sys_iter->second.MeanSeeing() << " " <<
	sys_iter->second.MeanExtinction() << " " <<
	sys_iter->second.MeanSky() << " " <<
	sys_iter->second.NObjects() << "\n";
    }
  }

  output_file.close();

  return 0;
}
