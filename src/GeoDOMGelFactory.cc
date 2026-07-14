#include <CLHEP/Units/SystemOfUnits.h>

#include <G4DisplacedSolid.hh>
#include <G4Orb.hh>
#include <G4RotationMatrix.hh>
#include <G4Sphere.hh>
#include <G4SubtractionSolid.hh>
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4VSolid.hh>
#include <GeoDOMGelFactory.hh>
#include <RAT/DB.hh>
#include <RAT/UnionSolidArray.hh>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace END {

namespace {
// Solve the 4x4 system M u = b in place (Gaussian elimination, partial pivot).
void Solve4(double M[4][4], double b[4], double u[4]) {
  for (int col = 0; col < 4; col++) {
    int piv = col;
    for (int r = col + 1; r < 4; r++)
      if (std::fabs(M[r][col]) > std::fabs(M[piv][col])) piv = r;
    if (piv != col) {
      for (int k = 0; k < 4; k++) std::swap(M[col][k], M[piv][k]);
      std::swap(b[col], b[piv]);
    }
    for (int r = 0; r < 4; r++) {
      if (r == col) continue;
      double f = M[r][col] / M[col][col];
      for (int k = col; k < 4; k++) M[r][k] -= f * M[col][k];
      b[r] -= f * b[col];
    }
  }
  for (int i = 0; i < 4; i++) u[i] = b[i] / M[i][i];
}
}  // namespace

// Algebraic sphere fit: minimise sum(|p-c|^2 - R^2)^2, which is linear in
// (cx, cy, cz, g=R^2-|c|^2):   [2x 2y 2z 1] . (cx,cy,cz,g) = x^2+y^2+z^2.
G4ThreeVector GeoDOMGelFactory::FitDOMCentre(const std::vector<double> &x, const std::vector<double> &y,
                                             const std::vector<double> &z, int lo, int hi) {
  double M[4][4] = {{0}}, rhs[4] = {0};
  for (int i = lo; i <= hi; i++) {
    double a[4] = {2 * x[i], 2 * y[i], 2 * z[i], 1.0};
    double bi = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
    for (int r = 0; r < 4; r++) {
      for (int c = 0; c < 4; c++) M[r][c] += a[r] * a[c];
      rhs[r] += a[r] * bi;
    }
  }
  double u[4];
  Solve4(M, rhs, u);
  return G4ThreeVector(u[0], u[1], u[2]);
}

G4VSolid *GeoDOMGelFactory::ConstructSolid(RAT::DBLinkPtr table) {
  std::string name = table->GetIndex();
  G4double r_min = table->GetD("r_min") * CLHEP::mm;  // photocathode-edge depth (~184)
  G4double r_max = table->GetD("r_max") * CLHEP::mm;  // inner glass (~203)
  int lo = table->GetI("start_idx");
  int hi = table->GetI("end_idx");

  // Dome radius + axial centre come from the PMT model itself: the widest ring
  // of the revolution profile is the hemisphere's equator, so its rho is the
  // bulb radius and its z is the dome centre. Reading them here keeps the gel in
  // lock-step with the PMT geometry (no duplicated numbers to drift).
  RAT::DBLinkPtr pmt = RAT::DB::Get()->GetLink("PMT", table->GetS("pmt_model"));
  const std::vector<double> &re = pmt->GetDArray("rho_edge");
  const std::vector<double> &ze = pmt->GetDArray("z_edge");
  int imax = std::max_element(re.begin(), re.end()) - re.begin();
  double dome_offset = ze[imax];                            // mm along +axis, PMT origin -> dome centre
  G4double r_dome = (re[imax] + 0.05) * CLHEP::mm;          // bulb radius + 0.05 mm gel/PMT clearance

  // Same PMTINFO table the pmtarray uses (global positions + radial directions).
  RAT::DBLinkPtr pos = RAT::DB::Get()->GetLink(table->GetS("pos_table"));
  const std::vector<double> &x = pos->GetDArray("x");
  const std::vector<double> &y = pos->GetDArray("y");
  const std::vector<double> &z = pos->GetDArray("z");
  const std::vector<double> &dx = pos->GetDArray("dir_x");
  const std::vector<double> &dy = pos->GetDArray("dir_y");
  const std::vector<double> &dz = pos->GetDArray("dir_z");

  // DOM centre = this volume's local origin (matches the dom_glass placement).
  G4ThreeVector centre = FitDOMCentre(x, y, z, lo, hi);

  G4VSolid *shell = new G4Sphere(name, r_min, r_max, 0., CLHEP::twopi, 0., CLHEP::pi);

  // A dome is a sphere -> isotropic -> only a translation is needed to carve it.
  std::vector<G4VSolid *> domes(hi - lo + 1);
  for (int i = lo; i <= hi; i++) {
    G4ThreeVector p(x[i], y[i], z[i]);
    G4ThreeVector dir = G4ThreeVector(dx[i], dy[i], dz[i]).unit();
    G4ThreeVector c = ((p - centre) + dome_offset * dir) * CLHEP::mm;  // dome centre, local frame
    domes[i - lo] = new G4DisplacedSolid(name + "_d" + std::to_string(i), new G4Orb(name + "_o", r_dome),
                                         G4Transform3D(G4RotationMatrix(), c));
  }

  return new G4SubtractionSolid(name + "_gel", shell, RAT::MakeUnionSolidArray(name + "_domes", domes));
}

}  // namespace END
