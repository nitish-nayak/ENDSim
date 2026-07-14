#ifndef __END_GeoDOMGelFactory__
#define __END_GeoDOMGelFactory__

#include <G4ThreeVector.hh>
#include <RAT/GeoSolidFactory.hh>
#include <vector>

namespace END {

// Optical-coupling gel for a KM3NeT DOM: a spherical shell (photocathode-edge
// depth -> inner glass) with the PMT domes subtracted -- the faithful meniscus
// that couples glass -> gel -> photocathode with no air gap.
//
// The JSON entry just names the same pos_table and pmt_model the pmtarray uses;
// the factory reads the PMT positions/directions and the PMT profile, and
// computes the dome radius, axial offset, and per-PMT dome centres itself:
//   r_min, r_max, start_idx, end_idx, pos_table, pmt_model
class GeoDOMGelFactory : public RAT::GeoSolidFactory {
 public:
  GeoDOMGelFactory() : GeoSolidFactory("domgel"){};
  virtual G4VSolid *ConstructSolid(RAT::DBLinkPtr table);

 private:
  // Least-squares sphere centre of PMT positions [lo, hi] = the DOM centre,
  // i.e. this volume's local-frame origin (matches the dom_glass placement).
  static G4ThreeVector FitDOMCentre(const std::vector<double> &x, const std::vector<double> &y,
                                    const std::vector<double> &z, int lo, int hi);
};

}  // namespace END

#endif
