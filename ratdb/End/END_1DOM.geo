{
name: "GEO",
index: "world",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "", // world volume has no mother
type: "box",
size: [250000.0, 250000.0, 250000.0], // mm, half-length
material: "water",
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
material: "lake_water",
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
material: "lake_water",
color: [0.0, 0.5, 1.0, 0.1],
invisible: 1,
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
pmt_model: "r12199",
mu_metal: 0,
mu_metal_material: "aluminum",
mu_metal_surface: "aluminum",
light_cone: 0,
pmt_detector_type: "idpmt",
sensitive_detector: "/mydet/pmt/inner",
efficiency_correction: 0.70000,
pos_table: "DOMINFO_aframe_spacing5m_hex_v3",
orientation: "manual",
invisible: 0,
color: [1.0, 0.0, 1.0, 0.5],
drawstyle: "wireframe"
}
