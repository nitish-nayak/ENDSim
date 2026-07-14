#ifndef __END_VertexGen_EcoMug__
#define __END_VertexGen_EcoMug__

#include <RAT/GLG4VertexGen.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

namespace END {

// Sea-level cosmic-muon vertex generator (ECoMUG, Pagano et al. NIM A 1014 (2021)).
// Samples (p, cos theta) from the horizontal-surface flux by Metropolis-Hastings
// (see EcoMugFlux.hh), azimuth uniform; emits one downward mu+/mu- per event.
// Vertical axis is y, downward is -y. The enclosing combo generator's position
// sets the injection-plane centre + height; this generator spreads each muon over
// a horizontal rectangle of half-size (fHalfX, fHalfZ) about that centre.
//
//   /generator/add combo ecomug:point:poisson
//   /generator/vtx/set [charge_ratio [half_x_mm [half_z_mm]]]
//
// At construction it logs the integrated area flux [Hz/m^2]; SetState logs the
// total rate for the configured plane, i.e. the value to hardcode into
// /generator/rate/set.
class VertexGen_EcoMug : public GLG4VertexGen {
 public:
  VertexGen_EcoMug(const char *arg_dbname = "ecomug");
  virtual ~VertexGen_EcoMug() {}
  virtual void GeneratePrimaryVertex(G4Event *event, G4ThreeVector &dx, G4double dt);
  virtual void SetState(G4String state);
  virtual G4String GetState();

 private:
  void Step();                    // advance the chain by NTHIN steps

  G4double fChargeRatio;          // N(mu+)/N(mu-)
  G4double fHalfX, fHalfZ;        // injection-plane half-size in x, z [mm]
  G4double fP, fCth, fPdf;        // persistent Metropolis chain state
};

}  // namespace END

#endif
