#include <VertexGen_EcoMug.hh>
#include <EcoMugFlux.hh>

#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <Randomize.hh>
#include <CLHEP/Units/SystemOfUnits.h>
#include <RAT/GLG4StringUtil.hh>
#include <RAT/Log.hh>

#include <sstream>

namespace END {

using namespace END::ecomug;

// Reach equilibrium before the first event (matches ecomug.f90 gen_horiz_init).
namespace {
constexpr int NWARMUP = 50000;
}

VertexGen_EcoMug::VertexGen_EcoMug(const char *arg_dbname)
    : GLG4VertexGen(arg_dbname), fChargeRatio(1.27), fHalfX(1500.0), fHalfZ(1500.0) {
  fP = 1.0;  fCth = 0.999;  fPdf = horpdf(fP, fCth);
  for (int i = 0; i < NWARMUP; ++i) Step();

  RAT::info << "VertexGen_EcoMug: integrated area flux = " << integrated_area_flux_hz()
            << " Hz/m^2 (p in [" << PMU_LO << ", " << PMU_HI << "] GeV/c)" << newline;
}

void VertexGen_EcoMug::Step() {
  mh_thinned(fP, fCth, fPdf, [] { return G4RandGauss::shoot(0.0, 1.0); }, [] { return G4UniformRand(); });
}

void VertexGen_EcoMug::GeneratePrimaryVertex(G4Event *event, G4ThreeVector &dx, G4double dt) {
  Step();

  double phi   = G4UniformRand() * CLHEP::twopi;
  double sinth = std::sqrt(std::max(0.0, 1.0 - fCth * fCth));
  double pmag  = fP * CLHEP::GeV;                      // GeV/c -> G4 units

  // Downward-going along -y; injection plane is the horizontal x-z plane.
  G4ThreeVector mom(pmag * sinth * std::cos(phi), -pmag * fCth, pmag * sinth * std::sin(phi));

  // dx (from the combo point pos-gen) is the plane centre; spread over the plane.
  G4ThreeVector pos = dx + G4ThreeVector((2.0 * G4UniformRand() - 1.0) * fHalfX, 0.0,
                                         (2.0 * G4UniformRand() - 1.0) * fHalfZ);

  int pdg = (G4UniformRand() < fChargeRatio / (1.0 + fChargeRatio)) ? -13 : 13;
  G4ParticleDefinition *mu = G4ParticleTable::GetParticleTable()->FindParticle(pdg);

  G4PrimaryVertex   *vertex   = new G4PrimaryVertex(pos, dt);
  G4PrimaryParticle *particle = new G4PrimaryParticle(mu, mom.x(), mom.y(), mom.z());
  particle->SetMass(MU_MASS_GEV * CLHEP::GeV);
  vertex->SetPrimary(particle);
  event->AddPrimaryVertex(vertex);
}

void VertexGen_EcoMug::SetState(G4String state) {
  state = util_strip_default(state);
  std::istringstream is(state);
  double cr, hx, hz;                        // "charge_ratio  half_x_mm  half_z_mm", all optional
  if (is >> cr) fChargeRatio = cr;
  if (is >> hx) fHalfX = hx;
  if (is >> hz) fHalfZ = hz;

  double area_m2 = (2.0 * fHalfX) * (2.0 * fHalfZ) / 1.0e6;   // mm^2 -> m^2
  RAT::info << "VertexGen_EcoMug: plane = " << area_m2 << " m^2  ->  set /generator/rate/set "
            << integrated_area_flux_hz() * area_m2 << " (Hz)" << newline;
}

G4String VertexGen_EcoMug::GetState() {
  std::ostringstream os;
  os << fChargeRatio << " " << fHalfX << " " << fHalfZ;
  return G4String(os.str());
}

}  // namespace END
