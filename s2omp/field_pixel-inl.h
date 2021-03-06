// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains a variant on the pixel class.  The goal here is to
// encode some manner of field field (galaxy density on the sky, CMB
// temperature, etc.) where we need both the value and an associated noise on
// the value.  Likewise, as will be seen in the field_union class, these pixels
// will generally be used to sample this field field uniformly over some
// region.  This is in contrast with the pixel_union object where the goal is to
// accurately describe a region's geometry using pixels of various sizes.

#ifndef SCALAR_PIXEL_H_
#define SCALAR_PIXEL_H_

#include "core.h"
#include "bound_interface.h"
#include "pixel.h"
#include "point.h"

namespace s2omp {

class field_pixel;

typedef std::vector<field_pixel> field_vector;
typedef field_vector::iterator field_iterator;
typedef field_vector::const_iterator field_const_iterator;
typedef std::pair<field_const_iterator, field_const_iterator> field_pair;
typedef std::vector<field_pixel *> field_ptr_vector;
typedef field_ptr_vector::const_iterator field_ptr_iterator;

class field_pixel: public pixel {
  // In order to do correlation function calculations, we need some
  // functionality beyond the normal pixel object.  In particular, we want
  // to be able to encode fields, which may take one of three forms:
  //
  // * Pure field quantities (e.g. CMB temperature or radio flux).
  // * Point-field densities (e.g. the projected galaxy density over some area).
  // * Point-sampled averages (e.g. the mean galaxy magnitude over some area).
  //
  // In order to accommodate those three cases, we need a field value (double)
  // and a point count (int), along with a weight (double).  pixels are taken
  // to be units of a pixel_union where the geometry of the map is the union of
  // all of the pixels and the pixels are not assumed to be at the same
  // resolution level.  field_pixels, OTOH, form the basis for a field_union
  // where the map is taken to be a regular sampling of some field over a given
  // area.  The total area for the map can be calculated, but operations like
  // determining whether or not a given position is inside or outside the
  // map is not generically available.
public:
  field_pixel(uint64 id);
  field_pixel(uint64 id, double intensity, double weight, long n_points);
  virtual ~field_pixel() {};

  inline double intensity() const {
    return intensity_;
  }
  inline double weight() const {
    return weight_;
  }
  inline long n_points() const {
    return n_points_;
  }

  inline void set_intensity(double intensity) {
    intensity_ = intensity;
  }
  inline void set_weight(double weight) {
    weight_ = weight;
  }
  inline void set_n_points(long n_points) {
    n_points_ = n_points;
  }

  inline double mean_intensity() const;
  inline void add_to_intensity(const double intensity, const long n_point);

  // Convert from our raw values to over-densities.  The return value indicates
  // whether or not the conversion was successful.
  inline bool convert_to_overdensity(double expected_intensity);
  inline bool convert_to_fractional_overdensity(double expected_intensity);

  // And two complementary methods to take us from over-densities back to raw
  // intensities.
  inline bool convert_from_overdensity(double expected_intensity);
  inline bool convert_from_fractional_overdensity(double expected_intensity);

  // Finally, a method to tell us whether we're in over-density mode or not.
  inline bool is_overdensity() const {
    return is_overdensity_;
  }

  // For the purposes of speeding up iterations in correlation functions, we
  // override the default behavior for get_center();
  inline virtual point get_center() const {
    return center_;
  }

private:
  field_pixel();

  double intensity_, weight_;
  point center_;
  long n_points_;
  bool is_overdensity_;
};

inline field_pixel::field_pixel() {
  set_id(0ULL);
  intensity_ = 0.0;
  weight_ = 0.0;
  n_points_ = 0L;
  center_ = point();
  is_overdensity_ = false;
}

inline field_pixel::field_pixel(uint64 id) {
  set_id(id);
  intensity_ = 0.0;
  weight_ = 0.0;
  n_points_ = 0L;
  center_ = center_point();
  is_overdensity_ = false;
}

inline field_pixel::field_pixel(
    uint64 id, double intensity, double weight, long n_points) {
  set_id(id);
  intensity_ = intensity;
  weight_ = weight;
  n_points_ = n_points;
  center_ = center_point();
  is_overdensity_ = false;
}

inline double field_pixel::mean_intensity() const {
  return n_points_ == 0 ? intensity_ : intensity_/n_points_;
}

inline void field_pixel::add_to_intensity(
    const double intensity, const long n_point) {
  intensity_ += intensity;
  n_points_ += n_point;
}

inline bool field_pixel::convert_to_overdensity(double expected_intensity) {
  if (is_overdensity_) {
    return false;
  }

  intensity_ -= expected_intensity;
  is_overdensity_ = true;
  return true;
}

inline bool field_pixel::convert_to_fractional_overdensity(
    double expected_intensity) {
  if (is_overdensity_) {
    return false;
  }

  intensity_ =
    (intensity_ - expected_intensity * weight_)/
    (expected_intensity * weight_);
  is_overdensity_ = true;
  return true;
}

inline bool field_pixel::convert_from_overdensity(double expected_intensity) {
  if (!is_overdensity_) {
    return false;
  }

  intensity_ = intensity_ + expected_intensity;
  is_overdensity_ = false;
  return true;
}

inline bool field_pixel::convert_from_fractional_overdensity(
    double expected_intensity) {
  if (!is_overdensity_) {
    return false;
  }

  double norm_intensity = expected_intensity * weight_;
  intensity_ = intensity_ * norm_intensity + norm_intensity;
  is_overdensity_ = false;
  return true;
}

} // end namespace s2omp

#endif /* SCALAR_PIXEL_H_ */
