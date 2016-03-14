// Copyright 2010  All Rights Reserved.
// Author: ryan.scranton@gmail.com (Ryan Scranton)

// STOMP is a set of libraries for doing astrostatistical analysis on the
// celestial sphere.  The goal is to enable descriptions of arbitrary regions
// on the sky which may or may not encode futher spatial information (galaxy
// density, CMB temperature, observational depth, etc.) and to do so in such
// a way as to make the analysis of that data as algorithmically efficient as
// possible.
//

#include "stomp_core.h"
#include "stomp_angular_bin.h"
#include "stomp_pixel.h"
#include "stomp_angular_coordinate.h"

namespace Stomp {

AngularBin::AngularBin() {
  theta_min_ = theta_max_ = sin2theta_min_ = sin2theta_max_ = 0.0;
  costheta_min_ = costheta_max_ = 1.0;
  weight_ = gal_gal_ = gal_rand_ = rand_gal_ = rand_rand_ = 0.0;
  pixel_wtheta_ = pixel_weight_ = wtheta_ = wtheta_error_ = 0.0;
  counter_ = 0;
  resolution_ = 0;
  ClearRegions();
  set_wtheta_ = false;
  set_wtheta_error_ = false;
}

AngularBin::AngularBin(double theta_min, double theta_max) {
  SetThetaMin(theta_min);
  SetThetaMax(theta_max);
  weight_ = gal_gal_ = gal_rand_ = rand_gal_ = rand_rand_ = 0.0;
  pixel_wtheta_ = pixel_weight_ = wtheta_ = wtheta_error_ = 0.0;
  counter_ = 0;
  ClearRegions();
  resolution_ = 0;
  set_wtheta_ = false;
  set_wtheta_error_ = false;
}

AngularBin::AngularBin(double theta_min, double theta_max, int16_t n_regions) {
  SetThetaMin(theta_min);
  SetThetaMax(theta_max);
  weight_ = gal_gal_ = gal_rand_ = rand_gal_ = rand_rand_ = 0.0;
  pixel_wtheta_ = pixel_weight_ = wtheta_ = wtheta_error_ = 0.0;
  counter_ = 0;
  ClearRegions();
  if (n_regions > 0) InitializeRegions(n_regions);
  resolution_ = 0;
  set_wtheta_ = false;
  set_wtheta_error_ = false;
}

AngularBin::~AngularBin() {
  theta_min_ = theta_max_ = sin2theta_min_ = sin2theta_max_ = 0.0;
  weight_ = gal_gal_ = gal_rand_ = rand_gal_ = rand_rand_ = 0.0;
  pixel_wtheta_ = pixel_weight_ = wtheta_ = wtheta_error_ = 0.0;
  counter_ = 0;
  ClearRegions();
  resolution_ = 0;
  set_wtheta_ = false;
  set_wtheta_error_ = false;
}

void AngularBin::ClearRegions() {
  weight_region_.clear();
  gal_gal_region_.clear();
  gal_rand_region_.clear();
  rand_gal_region_.clear();
  rand_rand_region_.clear();
  pixel_wtheta_region_.clear();
  pixel_weight_region_.clear();
  wtheta_region_.clear();
  wtheta_error_region_.clear();
  counter_region_.clear();
  n_region_ = 0;
}

void AngularBin::InitializeRegions(int16_t n_regions) {
  ClearRegions();
  if (n_regions > 0) {
  	n_region_ = n_regions;
  	n_region_pairs_ = n_region_ * n_region_;
    weight_region_.reserve(n_region_pairs_);
    gal_gal_region_.reserve(n_region_pairs_);
    gal_rand_region_.reserve(n_region_pairs_);
    rand_gal_region_.reserve(n_region_pairs_);
    rand_rand_region_.reserve(n_region_pairs_);
    pixel_wtheta_region_.reserve(n_region_pairs_);
    pixel_weight_region_.reserve(n_region_pairs_);
    wtheta_region_.reserve(n_region_pairs_);
    wtheta_error_region_.reserve(n_region_pairs_);
    counter_region_.reserve(n_region_pairs_);
    for (uint16_t k=0;k<n_region_pairs_;k++) {
      weight_region_[k] = 0.0;
      gal_gal_region_[k] = 0.0;
      gal_rand_region_[k] = 0.0;
      rand_gal_region_[k] = 0.0;
      rand_rand_region_[k] = 0.0;
      pixel_wtheta_region_[k] = 0.0;
      pixel_weight_region_[k] = 0.0;
      wtheta_region_[k] = 0.0;
      wtheta_error_region_[k] = 0.0;
      counter_region_[k] = 0;
  	}
  }
}

void AngularBin::SetResolution(uint32_t resolution) {
  resolution_ = resolution;
}

void AngularBin::CalculateResolution(double lammin, double lammax,
				     uint32_t max_resolution) {
  if (lammin < -70.0) {
    std::cout << "Stomp::AngularBin::CalculateResolution - " <<
      "Resetting minimum lambda value to -70.0...\n";
    lammin = -70.0;
  }
  if (lammax > 70.0) {
    std::cout << "Stomp::AngularBin::CalculateResolution - " <<
      "Resetting maximum lambda value to 70.0...\n";
    lammax = 70.0;
  }

  AngularCoordinate min_ang(lammin, 0.0, AngularCoordinate::Survey);
  AngularCoordinate max_ang(lammax, 0.0, AngularCoordinate::Survey);

  uint32_t pixel_resolution = HPixResolution;

  uint32_t ny_req = 1000000000;
  uint32_t small_good = 0, eta_good = 0;
  Pixel tmp_pix, tmp2_pix;

  while (((small_good < ny_req) || (eta_good < ny_req)) &&
	 (pixel_resolution <= max_resolution/2)) {

    small_good = eta_good = 0;
    pixel_resolution <<= 1;

    // Skip this resolution scale if the Lambda limits are smaller than pixels
    // at this resolution.
    if (lammax - lammin > sqrt(Pixel::Area(pixel_resolution))) {
      tmp_pix.SetResolution(pixel_resolution);
      tmp2_pix.SetResolution(pixel_resolution);

      tmp_pix.SetPixnumFromAng(min_ang);
      uint32_t ny_max = tmp_pix.PixelY();

      tmp_pix.SetPixnumFromAng(max_ang);
      uint32_t ny_min = tmp_pix.PixelY();

      ny_req = ny_max - ny_min;

      tmp2_pix = tmp_pix;
      for (uint32_t y=ny_min+1,x=tmp2_pix.PixelX();y<=ny_max;y++) {
        tmp_pix.SetPixnumFromXY(x+1,y);
        double costheta =
            tmp_pix.UnitSphereX()*tmp2_pix.UnitSphereX() +
            tmp_pix.UnitSphereY()*tmp2_pix.UnitSphereY() +
            tmp_pix.UnitSphereZ()*tmp2_pix.UnitSphereZ();
        if (1.0 - costheta*costheta < Sin2ThetaMax()) eta_good++;

        tmp_pix.SetPixnumFromXY(x,y);
        costheta =
            tmp_pix.UnitSphereX()*tmp2_pix.UnitSphereX() +
            tmp_pix.UnitSphereY()*tmp2_pix.UnitSphereY() +
            tmp_pix.UnitSphereZ()*tmp2_pix.UnitSphereZ();
        if (1.0 - costheta*costheta < Sin2ThetaMax()) small_good++;

        tmp2_pix = tmp_pix;
      }
    }
  }

  SetResolution(pixel_resolution);
}

void AngularBin::SetTheta(double theta) {
  theta_ = theta;
}

void AngularBin::SetThetaMin(double theta_min) {
  theta_min_ = theta_min;
  sin2theta_min_ =
    sin(theta_min_*DegToRad)*sin(theta_min_*DegToRad);
  costheta_max_ = cos(theta_min_*DegToRad);
}

void AngularBin::SetThetaMax(double theta_max) {
  theta_max_ = theta_max;
  sin2theta_max_ =
    sin(theta_max_*DegToRad)*sin(theta_max_*DegToRad);
  costheta_min_ = cos(theta_max_*DegToRad);
}

bool AngularBin::WithinBounds(double theta) {
  return (DoubleGE(theta, theta_min_) &&
	  DoubleLE(theta, theta_max_) ? true : false);
}

bool AngularBin::WithinSin2Bounds(double sin2theta) {
  return (DoubleGE(sin2theta, sin2theta_min_) &&
	  DoubleLE(sin2theta, sin2theta_max_) ? true : false);
}

bool AngularBin::WithinCosBounds(double costheta) {
  return (DoubleGE(costheta, costheta_min_) &&
	  DoubleLE(costheta, costheta_max_) ? true : false);
}

double AngularBin::Area() {
  return (costheta_max_ - costheta_min_)*2.0*Pi*StradToDeg;
}

double AngularBin::PoissonNoise(double objects_per_square_degree,
                                double survey_area) {
  return 1.0/sqrt(objects_per_square_degree*objects_per_square_degree*
                  survey_area*Area());
}

void AngularBin::AddToPixelWtheta(double dwtheta, double dweight,
				  int16_t region_a, int16_t region_b) {
  pixel_wtheta_ += dwtheta;
  pixel_weight_ += dweight;

  if ((region_a != -1) && (region_b != -1)) {
	  pixel_wtheta_region_[region_a * n_region_ + region_b] += dwtheta;
	  pixel_weight_region_[region_a * n_region_ + region_b] += dweight;
  }
}

void AngularBin::AddToWeight(double weight, int16_t region_a, int16_t region_b) {
  weight_ += weight;

  if ((region_a != -1) && (region_b != -1)) {
    weight_region_[region_a * n_region_ + region_b] += weight;
  }
}

void AngularBin::AddToCounter(uint32_t step, int16_t region_a, int16_t region_b) {
  counter_ += step;

  if ((region_a != -1) && (region_b != -1)) {
    counter_region_[region_a * n_region_ + region_b] += step;
  }
}

void AngularBin::MoveWeightToGalGal() {
  gal_gal_ += weight_;
  weight_ = 0.0;
  for (int16_t k=0;k<n_region_pairs_;k++) {
    gal_gal_region_[k] += weight_region_[k];
    weight_region_[k] = 0.0;
  }
}

void AngularBin::MoveWeightToGalRand(bool move_to_rand_gal) {
  gal_rand_ += weight_;
  if (move_to_rand_gal) rand_gal_ += weight_;
  weight_ = 0.0;
  for (int16_t k=0;k<n_region_pairs_;k++) {
    gal_rand_region_[k] += weight_region_[k];
    if (move_to_rand_gal) rand_gal_region_[k] += weight_region_[k];
    weight_region_[k] = 0.0;
  }
}

void AngularBin::MoveWeightToRandGal(bool move_to_gal_rand) {
  rand_gal_ += weight_;
  if (move_to_gal_rand) gal_rand_ += weight_;
  weight_ = 0.0;
  for (int16_t k=0;k<n_region_pairs_;k++) {
    rand_gal_region_[k] += weight_region_[k];
    if (move_to_gal_rand) gal_rand_region_[k] += weight_region_[k];
    weight_region_[k] = 0.0;
  }
}

void AngularBin::MoveWeightToRandRand() {
  rand_rand_ += weight_;
  weight_ = 0.0;
  for (int16_t k=0;k<n_region_pairs_;k++) {
    rand_rand_region_[k] += weight_region_[k];
    weight_region_[k] = 0.0;
  }
}

void AngularBin::RescaleGalGal(double weight) {
  gal_gal_ /= weight;
  for (int16_t k=0;k<n_region_pairs_;k++) gal_gal_region_[k] /= weight;
}

void AngularBin::RescaleGalRand(double weight) {
  gal_rand_ /= weight;
  for (int16_t k=0;k<n_region_pairs_;k++) gal_rand_region_[k] /= weight;
}

void AngularBin::RescaleRandGal(double weight) {
  rand_gal_ /= weight;
  for (int16_t k=0;k<n_region_pairs_;k++) rand_gal_region_[k] /= weight;
}

void AngularBin::RescaleRandRand(double weight) {
  rand_rand_ /= weight;
  for (int16_t k=0;k<n_region_pairs_;k++) rand_rand_region_[k] /= weight;
}

void AngularBin::Reset() {
  weight_ = gal_gal_ = gal_rand_ = rand_gal_ = rand_rand_ = 0.0;
  pixel_wtheta_ = pixel_weight_ = wtheta_ = wtheta_error_ = 0.0;
  counter_ = 0;
  if (n_region_ > 0) {
    for (int16_t k=0;k<n_region_pairs_;k++) {
      weight_region_[k] = 0.0;
      gal_gal_region_[k] = 0.0;
      gal_rand_region_[k] = 0.0;
      rand_gal_region_[k] = 0.0;
      rand_rand_region_[k] = 0.0;
      pixel_wtheta_region_[k] = 0.0;
      pixel_weight_region_[k] = 0.0;
      wtheta_region_[k] = 0.0;
      wtheta_error_region_[k] = 0.0;
      counter_region_[k] = 0;
    }
  }
}

void AngularBin::ResetPixelWtheta() {
  pixel_wtheta_ = 0.0;
  pixel_weight_ = 0.0;
  if (n_region_ > 0) {
    for (int16_t k=0;k<n_region_pairs_;k++) {
      pixel_wtheta_region_[k] = 0.0;
      pixel_weight_region_[k] = 0.0;
    }
  }
}

void AngularBin::ResetWeight() {
  weight_ = 0.0;
  if (n_region_ > 0)
    for (int16_t k=0;k<n_region_pairs_;k++) weight_region_[k] = 0.0;
}

void AngularBin::ResetCounter() {
  counter_ = 0;
  if (n_region_ > 0)
    for (int16_t k=0;k<n_region_pairs_;k++) counter_region_[k] = 0;
}

void AngularBin::ResetGalGal() {
  gal_gal_ = 0.0;
  if (n_region_ > 0)
    for (int16_t k=0;k<n_region_pairs_;k++) gal_gal_region_[k] = 0.0;
}

void AngularBin::ResetGalRand() {
  gal_rand_ = 0.0;
  if (n_region_ > 0)
    for (int16_t k=0;k<n_region_pairs_;k++) gal_rand_region_[k] = 0.0;
}

void AngularBin::ResetRandGal() {
  rand_gal_ = 0.0;
  if (n_region_ > 0)
    for (int16_t k=0;k<n_region_pairs_;k++) rand_gal_region_[k] = 0.0;
}

void AngularBin::ResetRandRand() {
  rand_rand_ = 0.0;
  if (n_region_ > 0)
    for (int16_t k=0;k<n_region_pairs_;k++) rand_rand_region_[k] = 0.0;
}

uint32_t AngularBin::Resolution() {
  return resolution_;
}

int16_t AngularBin::NRegion() {
  return n_region_;
}

double AngularBin::Theta() {
  return theta_;
}

double AngularBin::ThetaMin() {
  return theta_min_;
}

double AngularBin::ThetaMax() {
  return theta_max_;
}

double AngularBin::Sin2ThetaMin() {
  return sin2theta_min_;
}

double AngularBin::Sin2ThetaMax() {
  return sin2theta_max_;
}

double AngularBin::CosThetaMin() {
  return costheta_min_;
}

double AngularBin::CosThetaMax() {
  return costheta_max_;
}

double AngularBin::Wtheta(int16_t region) {
  if (set_wtheta_) {
    return (region == -1 ? wtheta_ :
	    (region < n_region_pairs_ ? wtheta_region_[region] : -1.0));
  } else {
    if (resolution_ == 0) {
    	if (region == -1) {
    		return (gal_gal_ - gal_rand_ - rand_gal_ + rand_rand_)/rand_rand_;
    	} else {
    		double total = 0;
    		double rand_rand = 0;
    		for (int16_t k=0;k<n_region_;k++) {
    			for (int16_t l=0;l<n_region_;l++) {
    				if ((k != region) && (l != region)) {
    					total += (gal_gal_region_[k * n_region_ + l] -
    							      gal_rand_region_[k * n_region_ + l] -
    							      rand_gal_region_[k * n_region_ + l] +
    							      rand_rand_region_[k * n_region_ + l]);
    					rand_rand += rand_rand_region_[k * n_region_ + l];
    				}
    			}
    		}
    		return total / (1.0 * rand_rand);
    	}
    } else {
    	if (region == -1) {
        return pixel_wtheta_/pixel_weight_;
    	} else {
    		double tot_wtheta = 0;
    		double tot_weight = 0;
    		for (int16_t k=0;k<n_region_;k++) {
    		  for (int16_t l=0;l<n_region_;l++) {
    		  	if ((k != region) && (l != region)) {
    		  	  tot_wtheta += pixel_wtheta_region_[k * n_region_ + l];
    		  	  tot_weight += pixel_weight_region_[k * n_region_ + l];
    		  	}
    		  }
    		}
    		return tot_wtheta / tot_weight;
    	}
    }
  }
}

double AngularBin::WthetaError(int16_t region) {
  if (set_wtheta_error_) {
    return (region == -1 ? wtheta_error_ :
	    (region < n_region_pairs_ ? wtheta_error_region_[region] : -1.0));
  } else {
    if (resolution_ == 0) {
      return (region == -1 ? 1.0/sqrt(gal_gal_) :
	      (region < n_region_pairs_ ?
	       1.0/sqrt(gal_gal_region_[region]) : -1.0));
    } else {
      return (region == -1 ? 1.0/sqrt(pixel_weight_) :
	      (region < n_region_pairs_ ? 1.0/sqrt(pixel_weight_region_[region]) :
	       -1.0));
    }
  }
}

double AngularBin::WeightedCrossCorrelation(int16_t region) {
  return (region == -1 ? weight_/counter_ :
	  (region < n_region_pairs_ ?
	   weight_region_[region]/counter_region_[region] : -1.0));
}

double AngularBin::PixelWtheta(int16_t region) {
  return (region == -1 ? pixel_wtheta_ :
	  (region < n_region_pairs_ ? pixel_wtheta_region_[region] : -1.0));
}

double AngularBin::PixelWeight(int16_t region) {
  return (region == -1 ? pixel_weight_ :
	  (region < n_region_pairs_ ? pixel_weight_region_[region] : -1.0));
}

double AngularBin::Weight(int16_t region) {
  return (region == -1 ? weight_ :
	  (region < n_region_pairs_ ? weight_region_[region] : -1.0));
}

uint32_t AngularBin::Counter(int16_t region) {
  return (region == -1 ? counter_ :
	  (region < n_region_pairs_ ? counter_region_[region] : -1));
}

double AngularBin::GalGal(int16_t region) {
  return (region == -1 ? gal_gal_ :
	  (region < n_region_pairs_ ? gal_gal_region_[region] : -1.0));
}

double AngularBin::GalRand(int16_t region) {
  return (region == -1 ? gal_rand_ :
	  (region < n_region_pairs_ ? gal_rand_region_[region] : -1.0));
}

double AngularBin::RandGal(int16_t region) {
  return (region == -1 ? rand_gal_ :
	  (region < n_region_pairs_ ? rand_gal_region_[region] : -1.0));
}

double AngularBin::RandRand(int16_t region) {
  return (region == -1 ? rand_rand_ :
	  (region < n_region_pairs_ ? rand_rand_region_[region] : -1.0));
}

double AngularBin::MeanWtheta() {
  double mean_wtheta = 0.0;
  for (int16_t k=0;k<n_region_;k++)
    mean_wtheta += Wtheta(k)/(1.0*n_region_);
  return mean_wtheta;
}

double AngularBin::MeanWthetaError() {
  double mean_wtheta = MeanWtheta();
  double mean_wtheta_error = 0.0;
  for (int16_t k=0;k<n_region_;k++)
    mean_wtheta_error += (mean_wtheta - Wtheta(k)) *
      (mean_wtheta - Wtheta(k));
  return (n_region_ == 0 ? 0.0 :
	  (n_region_ - 1.0) / (1.0 * n_region_) * sqrt(mean_wtheta_error));
}

double AngularBin::MeanWeightedCrossCorrelation() {
  double mean_weight_cross_correlation = 0.0;
  for (int16_t k=0;k<n_region_;k++)
    mean_weight_cross_correlation +=
      1.0*weight_region_[k]/counter_region_[k]/(1.0*n_region_);
  return mean_weight_cross_correlation;
}

double AngularBin::MeanWeightedCrossCorrelationError() {
  double mean_weighted_cross_correlation = MeanWeightedCrossCorrelation();
  double mean_weighted_cross_correlation_error = 0.0;
  for (int16_t k=0;k<n_region_;k++)
    mean_weighted_cross_correlation_error +=
      (mean_weighted_cross_correlation - WeightedCrossCorrelation(k))*
      (mean_weighted_cross_correlation - WeightedCrossCorrelation(k));
  return (n_region_ == 0 ? 0.0 :
	  sqrt(mean_weighted_cross_correlation_error)/(1.0 * n_region_));
}

double AngularBin::MeanWeight() {
  double mean_weight = 0.0;
  for (int16_t k=0;k<n_region_;k++)
    mean_weight += weight_region_[k]/(1.0*n_region_);
  return mean_weight;
}

double AngularBin::MeanCounter() {
  double mean_counter = 0.0;
  for (int16_t k=0;k<n_region_;k++)
    mean_counter += 1.0*counter_region_[k]/(1.0*n_region_);
  return mean_counter;
}

double AngularBin::MeanGalGal() {
  double mean_gal_gal = 0.0;
  for (int16_t k=0;k<n_region_;k++)
    mean_gal_gal += gal_gal_region_[k]/(1.0*n_region_);
  return mean_gal_gal;
}

double AngularBin::MeanGalRand() {
  double mean_gal_rand = 0.0;
  for (int16_t k=0;k<n_region_;k++)
    mean_gal_rand += gal_rand_region_[k]/(1.0*n_region_);
  return mean_gal_rand;
}

double AngularBin::MeanRandGal() {
  double mean_rand_gal = 0.0;
  for (int16_t k=0;k<n_region_;k++)
    mean_rand_gal += rand_gal_region_[k]/(1.0*n_region_);
  return mean_rand_gal;
}

double AngularBin::MeanRandRand() {
  double mean_rand_rand = 0.0;
  for (int16_t k=0;k<n_region_;k++)
    mean_rand_rand += rand_rand_region_[k]/(1.0*n_region_);
  return mean_rand_rand;
}

bool AngularBin::ThetaOrder(AngularBin theta_a, AngularBin theta_b) {
  return (theta_a.ThetaMin() < theta_b.ThetaMin() ? true : false);
}

bool AngularBin::SinThetaOrder(AngularBin theta_a, AngularBin theta_b) {
  return (theta_a.Sin2ThetaMin() < theta_b.Sin2ThetaMin() ? true : false);
}

bool AngularBin::ReverseResolutionOrder(AngularBin theta_a,
					AngularBin theta_b) {
  return (theta_b.Resolution() < theta_a.Resolution() ? true : false);
}

} // end namespace Stomp
