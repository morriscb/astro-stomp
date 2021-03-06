// Copyright 2012  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//
// This header file contains the pixel_union class.  A pixel_union is analogous
// to the S2::CellUnion or Stomp::Map classes, where we describe a region of
// the sky with a collection of pixels.  Most of the functionality of the
// Map class is replicated here, with the exception of the fact that
// pixel_union objects do not encode scalar weights.

#ifndef PIXEL_UNION_H_
#define PIXEL_UNION_H_

#include <stdint.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include "core.h"
#include "bound_interface.h"
#include "circle_bound.h"
#include "pixel.h"

namespace s2omp {

class point;

class pixel_union : public bound_interface {
public:
  // Like the S2CellUnion and Stomp::Map classes, the pixels that comprise
  // a pixel_union need to be properly normalized to make the remainder of
  // the class work.  Since that can be time-consuming, we defer that
  // operation to the init() function.
  pixel_union();
  virtual ~pixel_union();

  // Alternatively, one can use this static method to return a pointer
  // to a normalized pixel_union from a vector of pixels.
  inline static pixel_union* from_covering(pixel_vector& pixels);

  // In order to efficiently store and search the area covered by a given set
  // of pixels, we need to sort them by index and combine child pixels into
  // parent pixels where possible.  normalize() performs this task and we make
  // it a static method so that this capability is available without
  // initializing a full pixel_union.
  static void normalize(pixel_vector* pixels);

  // Use the input vector of pixels to initialize this pixel_union.  These
  // pixels will be normalized and stored internally.  The input pixel_vector
  // will be emptied by this method.
  void init(pixel_vector& pixels);

  // In some cases, we may wish to soften the edges of our pixel_union (for
  // instance if the current union is formed from combining multiple high
  // resolution unions producing a union that is unwieldy).  In this case,
  // we combine all pixels smaller than the input level, producing a parent
  // pixel at max_level if more than half its area was in the original union.
  void soften(int max_level);

  // As with Stomp::Map, we have multiple methods for doing logical operations
  // on pixel_unions: union, intersection and exclusion.  In all cases, the
  // contents of the current pixel_union are replaced by the results of the
  // operation, which may lead to an empty pixel_union.
  void combine_with(const pixel_union& pix_union);
  void intersect_with(const pixel_union& pix_union);
  void exclude_from(const pixel_union& pix_union);

  // Alternatively, we can initialize this pixel_union from the results of
  // doing a union, intersection or exclusion between two other pixel_unions.
  void init_from_combination(const pixel_union& a, const pixel_union& b);
  void init_from_intersection(const pixel_union& a, const pixel_union& b);
  void init_from_exclusion(const pixel_union& a, const pixel_union& b);

  bool intersects(const pixel& pix) const;
  bool may_intersect(const pixel_union& pix_union) const;
  bool intersects(const pixel_union& pix_union) const;

  // Method for returning the child pixels that intersect with a union
  void pixel_intersection(const pixel& pix, pixel_vector* pixels) const;
  static void pixel_intersection(const pixel_union& pix_union, const pixel& pix,
                                 pixel_vector* pixels);

  // Same for returning the parts of the input pixel that are outside a union.
  void pixel_exclusion(const pixel& pix, pixel_vector* pixels) const;
  static void pixel_exclusion(const pixel_union& pix_union, const pixel& pix,
                       pixel_vector* pixels);

  inline int min_level() const {
    return min_level_;
  }
  inline int max_level() const {
    return max_level_;
  }

  // Like with pixel, we can define a range_min and range_max of the pixels
  // that bound our pixel_union.  For compact regions, this can be used for
  // quick and dirty random point generation.
  inline pixel range_min() const {
    return range_min_;
  }
  inline pixel range_max() const {
    return range_max_;
  }

  inline pixel_iterator begin() const {
    return pixels_.begin();
  }
  inline pixel_iterator end() const {
    return pixels_.end();
  }

  // API from bound_interface.h
  inline virtual bool is_empty() const {
    return pixels_.empty();
  }
  inline virtual long size() const {
    return pixels_.size();
  }
  virtual void clear();
  inline virtual double area() const {
    return area_;
  }

  virtual bool contains(const point& p) const;
  virtual bool contains(const pixel& pix) const;

  virtual double contained_area(const pixel& pix) const;
  virtual bool may_intersect(const pixel& pix) const;

  virtual circle_bound get_bound() const;
  virtual point get_center() const;

  virtual void get_covering(pixel_vector* pixels) const;
  virtual void get_size_covering(
      const long max_pixels, pixel_vector* pixels) const;
  virtual void get_area_covering(
      double fractional_area_tolerance, pixel_vector* pixels) const;
  virtual void get_interior_covering(int max_level, pixel_vector* pixels) const;
  virtual void get_simple_covering(int level, pixel_vector* pixels) const;
  virtual void get_center_covering(int level, pixel_vector* pixels) const;


private:
  void generate_basic_covering(int level);
  void initialize_bound();

  // For automatic regionation is is useful to know what the average area of
  // a pixel is in our union.

  pixel_vector pixels_;
  circle_bound bound_;
  pixel range_min_;
  pixel range_max_;
  int min_level_, max_level_;
  double area_;
  bool initialized_, initialized_bound_;
};

inline pixel_union* pixel_union::from_covering(pixel_vector& pixels) {
  pixel_union* pix_union = new pixel_union();
  pix_union->init(pixels);
  return pix_union;
}

} // end namespace s2omp

#endif /* PIXEL_UNION_H_ */
