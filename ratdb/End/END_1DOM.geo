{
name: "GEO",
index: "world",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "", // world volume has no mother
type: "box",
size: [250000.0, 250000.0, 250000.0], // mm, half-length
material: "air", // surface setup: 1 DOM in a dark box
invisible: 0,
color: [1.0, 1.0, 1.0, 0.1],
drawstyle: "wireframe"
}

////////////////////////////////// Define detector properties.  ///////////////////////////////////
{
name: "GEO",
index: "detector",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "world",
type: "box",
size: [25000.0, 25000.0, 250000.0], // mm, half-length (50m x 50m x 500m total)
position: [0.0, 0.0, 0.0],
material: "air", // dark-box air (was lake_water)
invisible: 0,
color: [0.0, 0.5, 1.0, 0.2],
drawstyle: "wireframe"
}

{
name: "GEO",
index: "dom_glass_3",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "detector",
type: "sphere",
r_max: 216.0,
position: [-1500.0000, 0.0001, 2001.9238],
material: "borosilicate_glass",
color: [0.8, 0.9, 1.0, 0.3],
drawstyle: "wireframe"
}

{
name: "GEO",
index: "dom_inner_3",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "dom_glass_3",
type: "sphere",
r_max: 203.0,
position: [0.0, 0.0, 0.0],
material: "air", // DOM interior: air (optically dead); the gel handles glass -> PMT coupling
color: [0.0, 0.5, 1.0, 0.1],
invisible: 1,
drawstyle: "wireframe"
}

////////////////////////////////// Optical-coupling gel meniscus.  ///////////////////////////////////
// Shell from the photocathode-edge depth to the inner glass, with the PMT domes
// subtracted (built by the domgel factory from the PMT positions + model). Couples
// glass -> gel -> photocathode with no air gap.
{
name: "GEO",
index: "dom_gel_3",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "dom_inner_3",
type: "domgel",
position: [0.0, 0.0, 0.0],
r_max: 203.0, // mm, inner glass surface (r_min is computed by the factory from the Ø72 cathode edge)
pos_table: "DOMINFO_aframe_spacing5m_hex_latest",
pmt_model: "r12199_rev",
start_idx: 93,
end_idx: 123,
material: "borosilicate_glass",
color: [1.0, 1.0, 0.6, 0.35],
drawstyle: "wireframe"
}

{
name: "GEO",
index: "dom_pmts_3",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "dom_inner_3",
type: "pmtarray",
start_idx: 93,
end_idx: 123,
pmt_model: "r12199_rev",
mu_metal: 0,
mu_metal_material: "aluminum",
mu_metal_surface: "aluminum",
light_cone: 0,
pmt_detector_type: "idpmt",
sensitive_detector: "/mydet/pmt/inner",
efficiency_correction: 0.70000,
pos_table: "DOMINFO_aframe_spacing5m_hex_latest",
orientation: "manual",
invisible: 0,
color: [1.0, 0.0, 1.0, 0.5],
drawstyle: "wireframe"
}
