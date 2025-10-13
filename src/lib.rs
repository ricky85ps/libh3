extern crate static_assertions;

#[cfg(not(docsrs))]
mod libh3_sys {
    #![allow(improper_ctypes)]
    #![allow(unused)]
    #![allow(non_upper_case_globals)]
    #![allow(non_camel_case_types)]
    #![allow(non_snake_case)]
    include!("libh3_sys.rs");
}
use std::convert::From;
use std::mem::MaybeUninit;

#[cfg(docsrs)]
static_assertions::const_assert_eq!(libh3_sys::H3_VERSION_MAJOR, 4);
#[cfg(docsrs)]
static_assertions::const_assert!(libh3_sys::H3_VERSION_MINOR >= 2);

#[derive(Debug)]
pub enum H3Err {
    Success = 0,      // Success (no error)
    Failed,           // The operation failed but a more specific error is not available
    Domain, // Argument was outside of acceptable range (when a more specific error code is not available)
    LatLngDomain, // Latitude or longitude arguments were outside of acceptable range
    ResDomain, // Resolution argument was outside of acceptable range
    CellInvalid, // `H3Index` cell argument was not valid
    DirEdgeInvalid, // `H3Index` directed edge argument was not valid
    UndirEdgeInvalid, // `H3Index` undirected edge argument was not valid
    VertexInvalid, // `H3Index` vertex argument was not valid
    Pentagon, // Pentagon distortion was encountered which the algorithm could not handle it
    DuplicateInput, // Duplicate input was encountered in the arguments and the algorithm could not handle it
    NotNeighbors,   // `H3Index` cell arguments were not neighbors
    ResMismatch,    // `H3Index` cell arguments had incompatible resolutions
    MemoryAlloc,    // Necessary memory allocation failed
    MemoryBounds,   // Bounds of provided memory were not large enough
    OptionInvalid,  // Mode or flags argument was not valid.
    InvalidErrCode, // H3ErrorCodes from libh3 is invalid
}

#[allow(non_snake_case, unused_variables, unreachable_patterns)]
impl From<libh3_sys::H3Error> for H3Err {
    fn from(value: libh3_sys::H3Error) -> Self {
        match value {
            H3ErrorCodes_E_SUCCESS => H3Err::Success,
            H3ErrorCodes_E_FAILED => H3Err::Failed,
            H3ErrorCodes_E_DOMAIN => H3Err::Domain,
            H3ErrorCodes_E_LATLNG_DOMAIN => H3Err::LatLngDomain,
            H3ErrorCodes_E_RES_DOMAIN => H3Err::ResDomain,
            H3ErrorCodes_E_CELL_INVALID => H3Err::CellInvalid,
            H3ErrorCodes_E_DIR_EDGE_INVALID => H3Err::DirEdgeInvalid,
            H3ErrorCodes_E_UNDIR_EDGE_INVALID => H3Err::UndirEdgeInvalid,
            H3ErrorCodes_E_VERTEX_INVALID => H3Err::VertexInvalid,
            H3ErrorCodes_E_PENTAGON => H3Err::Pentagon,
            H3ErrorCodes_E_DUPLICATE_INPUT => H3Err::DuplicateInput,
            H3ErrorCodes_E_NOT_NEIGHBORS => H3Err::NotNeighbors,
            H3ErrorCodes_E_RES_MISMATCH => H3Err::ResMismatch,
            H3ErrorCodes_E_MEMORY_ALLOC => H3Err::MemoryAlloc,
            H3ErrorCodes_E_MEMORY_BOUNDS => H3Err::MemoryBounds,
            H3ErrorCodes_E_OPTION_INVALID => H3Err::OptionInvalid,
            _else => H3Err::InvalidErrCode,
        }
    }
}

/// Convert degrees to radians
///
/// ```
/// use libh3::degs_to_rads;
/// assert_eq!(2.413790355508158, degs_to_rads(138.3));
/// ```
pub fn degs_to_rads(degrees: f64) -> f64 {
    unsafe { libh3_sys::degsToRads(degrees) }
}

/// Convert radians to degrees
///
/// ```
/// use libh3::rads_to_degs;
/// assert_eq!(138.3, rads_to_degs(2.413790355508158));
/// ```
pub fn rads_to_degs(radians: f64) -> f64 {
    unsafe { libh3_sys::radsToDegs(radians) }
}

/// Represent a coordinate
#[derive(Debug, Clone)]
pub struct GeoCoord {
    // The latitute of the coordinate, typcially this should be specified using
    // radians but it is easy to convert using [degs_to_rads](degs_to_rads)
    pub lat: f64,
    // The longgitude of the coordinate, typcially this should be specified using
    // radians but it is easy to convert using [degs_to_rads](degs_to_rads)
    pub lng: f64,
}

impl GeoCoord {
    /// Create a new GeoCoord representing a coordinate
    ///
    /// # Arguments
    ///
    /// * `lat` - The latitude of the coordinate
    /// * `lng` - The longgitude of the coordinate
    ///
    pub fn new(lat: f64, lng: f64) -> GeoCoord {
        GeoCoord { lat, lng }
    }
}

impl From<libh3_sys::LatLng> for GeoCoord {
    fn from(coord: libh3_sys::LatLng) -> Self {
        GeoCoord {
            lat: coord.lat,
            lng: coord.lng,
        }
    }
}

impl From<&GeoCoord> for libh3_sys::LatLng {
    fn from(coord: &GeoCoord) -> Self {
        libh3_sys::LatLng {
            lat: coord.lat,
            lng: coord.lng,
        }
    }
}

/// A H3 index value a unique address of a hexagon or more unlikely
/// a pentagon.
pub type H3Index = libh3_sys::H3Index;

/// A resolution that ranges from 0 to 15.
///
/// See the [resolution table](https://h3geo.org/docs/core-library/restable)
/// for the sizes of resolution.
pub type Resolution = u8;

/// Return the edge length of a hexagon at a particular resolution in kilometers.
///
/// ```
/// use libh3::edge_length_km;
/// assert_eq!(edge_length_km(5)?, 9.85409099);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn edge_length_km(resolution: Resolution) -> Result<f64, H3Err> {
    let mut out = 0f64;
    match unsafe {
        libh3_sys::getHexagonEdgeLengthAvgKm(resolution as ::std::os::raw::c_int, &mut out)
    }
    .into()
    {
        H3Err::Success => Ok(out),
        error => Err(error),
    }
}

/// Return the number of hexagons at a particular resolution.
///
/// ```
/// use libh3::num_cells;
/// assert_eq!(num_cells(5)?, 2016842);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn num_cells(resolution: Resolution) -> Result<i64, H3Err> {
    let mut out = 0i64;
    match unsafe { libh3_sys::getNumCells(resolution as ::std::os::raw::c_int, &mut out) }.into() {
        H3Err::Success => Ok(out),
        error => Err(error),
    }
}

/// Return the edge length of a hexagon at a particular resolution in meters.
///
/// ```
/// use libh3::edge_length_m;
/// assert_eq!(edge_length_m(5)?, 9854.09099);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn edge_length_m(resolution: Resolution) -> Result<f64, H3Err> {
    let mut out = 0f64;
    match unsafe {
        libh3_sys::getHexagonEdgeLengthAvgM(resolution as ::std::os::raw::c_int, &mut out)
    }
    .into()
    {
        H3Err::Success => Ok(out),
        error => Err(error),
    }
}

/// Convert a GeoCoord to a H3 index.
///
/// ```
/// use libh3::{GeoCoord, degs_to_rads, geo_to_h3};
/// let coords = GeoCoord {
///   lat: degs_to_rads(40.689167),
///   lng: degs_to_rads(-74.044444),
/// };
///
/// let v = geo_to_h3(&coords, 10);
/// assert_eq!(v?, 0x8a2a1072b59ffff);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn geo_to_h3(coord: &GeoCoord, resolution: Resolution) -> Result<H3Index, H3Err> {
    let mut out = 0u64;
    match unsafe {
        libh3_sys::latLngToCell(&coord.into(), resolution as ::std::os::raw::c_int, &mut out)
    }
    .into()
    {
        H3Err::Success => Ok(out),
        error => Err(error),
    }
}

/// Convert a H3 index value to a GeoCoord
///
/// ```
/// use libh3::h3_to_geo;
/// let r = h3_to_geo(0x8a2a1072b59ffff)?;
/// assert_eq!(r.lat, 0.710164381905454);
/// assert_eq!(r.lng, -1.2923191206954798);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn h3_to_geo(h3: H3Index) -> Result<GeoCoord, H3Err> {
    let mut out = libh3_sys::LatLng {
        lat: 0f64,
        lng: 0f64,
    };
    match unsafe { libh3_sys::cellToLatLng(h3, &mut out) }.into() {
        H3Err::Success => Ok(out.into()),
        error => Err(error),
    }
}

pub struct CellBoundary {
    pub vec: Vec<GeoCoord>,
}

impl From<libh3_sys::CellBoundary> for CellBoundary {
    fn from(value: libh3_sys::CellBoundary) -> Self {
        let mut cellboundary = CellBoundary { vec: Vec::new() };
        for i in 0..value.numVerts as usize {
            cellboundary.vec.push(value.verts[i].into());
        }
        cellboundary
    }
}

/// Convert a H3 index value to a GeoBoundary which are a
/// vector of points that describe a H3 Index's boundary
///
/// ```
/// use libh3::h3_to_geo_boundary;
/// let foo = h3_to_geo_boundary(0x8a2a1072b59ffff)?;
/// assert_eq!(foo.vec.len(), 6);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn h3_to_geo_boundary(h3: H3Index) -> Result<CellBoundary, H3Err> {
    let mut boundary_result: MaybeUninit<libh3_sys::CellBoundary> =
        unsafe { MaybeUninit::uninit().assume_init() };

    match unsafe { libh3_sys::cellToBoundary(h3, boundary_result.as_mut_ptr()) }.into() {
        H3Err::Success => {}
        error => return Err(error),
    }
    let boundary = unsafe { boundary_result.assume_init() };
    let result = boundary.verts[0..boundary.numVerts as usize]
        .iter()
        .map(|&vert| GeoCoord::new(vert.lat, vert.lng))
        .collect();
    Ok(CellBoundary { vec: result })
}

/// Return the resolution of a H3 index
///
/// ```
/// use libh3::{GeoCoord, degs_to_rads, geo_to_h3, h3_get_resolution };
/// let coords = GeoCoord::new(
///   degs_to_rads(40.689167),
///   degs_to_rads(-74.044444),
/// );
///
/// let v = geo_to_h3(&coords, 10);
/// assert_eq!(h3_get_resolution(v.unwrap()), 10);
/// # Ok::<(), libh3::H3Err>(())
/// ````
pub fn h3_get_resolution(h3: H3Index) -> Resolution {
    unsafe { libh3_sys::getResolution(h3) as Resolution }
}

/// Determine if H3 index is valid
///
/// ```
/// use libh3::{GeoCoord, degs_to_rads, geo_to_h3, h3_is_valid};
/// let coords = GeoCoord {
///   lat: degs_to_rads(40.689167),
///   lng: degs_to_rads(-74.044444),
/// };
///
/// let v = geo_to_h3(&coords, 10);
/// assert_eq!(h3_is_valid(v.unwrap()),true);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn h3_is_valid(h3: H3Index) -> bool {
    (unsafe { libh3_sys::isValidCell(h3) } != 0)
}

/// Determine if two H3 indexes are neighbors
///
/// ```
/// use libh3::{GeoCoord, degs_to_rads, geo_to_h3, h3_indexes_are_neighbors};
/// let coords = GeoCoord {
///   lat: degs_to_rads(40.689167),
///   lng: degs_to_rads(-74.044444),
/// };
///
/// let v = geo_to_h3(&coords, 10)?;
/// assert_eq!(h3_indexes_are_neighbors(v, v)?, false);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn h3_indexes_are_neighbors(origin: H3Index, destination: H3Index) -> Result<bool, H3Err> {
    let mut out: ::std::os::raw::c_int = 0;
    match unsafe { libh3_sys::areNeighborCells(origin, destination, &mut out) }.into() {
        H3Err::Success => Ok(out != 0),
        error => Err(error),
    }
}

/// Determine the area of a hexagon at a particular resolution.
///
/// ```
/// use libh3::{ hex_area_km_2};
/// assert_eq!(hex_area_km_2(10)?, 0.01504750190766435);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn hex_area_km_2(resolution: i32) -> Result<f64, H3Err> {
    let mut out = 0f64;
    match unsafe { libh3_sys::getHexagonAreaAvgKm2(resolution as ::std::os::raw::c_int, &mut out) }
        .into()
    {
        H3Err::Success => Ok(out),
        error => Err(error),
    }
}

/// Determine if the specified H3 index is a pentagon.
/// ```
/// assert_eq!(libh3::h3_is_pentagon(0x8a2a1072b59ffff), false);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn h3_is_pentagon(h3: H3Index) -> bool {
    (unsafe { libh3_sys::isPentagon(h3) } != 0)
}

/// Get the number of the base cell for a given H3 index
///
/// ```
/// use libh3;
/// assert_eq!(libh3::h3_get_base_cell(0x8a2a1072b59ffff), 21);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn h3_get_base_cell(h3: H3Index) -> i32 {
    unsafe { libh3_sys::getBaseCellNumber(h3) }
}

pub type HexDistance = i32;

/// Get all hexagons in a k-ring around a given center. The order of the hexagons is undefined.
///
/// # Arguments
///
/// * `origin` - The center of the ring.
/// * `radius` - The radis of the ring in hexagons, which is the same resolution as the origin.
///
/// ```
/// let expected_kring = vec![
///   0x8a2a1072b59ffff,
///   0x8a2a1072b597fff,
///   0x8a2a1070c96ffff,
///   0x8a2a1072b4b7fff,
///   0x8a2a1072b4a7fff,
///   0x8a2a1072b58ffff,
///   0x8a2a1072b587fff,
/// ];
/// let r = libh3::k_ring(0x8a2a1072b59ffff, 1)?;
/// assert_eq!(r, expected_kring);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn k_ring(origin: H3Index, radius: HexDistance) -> Result<Vec<H3Index>, H3Err> {
    let mut max = 0;
    match unsafe { libh3_sys::maxGridDiskSize(radius, &mut max) }.into() {
        H3Err::Success => {}
        error => return Err(error),
    }

    let mut r: Vec<H3Index> = vec![0; max as usize];
    match unsafe { libh3_sys::gridDisk(origin, radius, r.as_mut_ptr()) }.into() {
        H3Err::Success => {}
        error => return Err(error),
    }
    r.retain(|v| *v != 0);
    Ok(r)
}

/// Get all hexagons in a k-ring around a given center, in an array of arrays
/// ordered by distance from the origin. The order of the hexagons within each ring is undefined.
///
/// # Arguments
///
/// * `origin` - The center of the ring.
/// * `radius` - The radis of the ring in hexagons, which is the same resolution as the origin.
///
/// ```
/// let expected_kring_distances = vec![
///   (0x8a2a1072b59ffff, 0),
///   (0x8a2a1072b597fff, 1),
///   (0x8a2a1070c96ffff, 1),
///   (0x8a2a1072b4b7fff, 1),
///   (0x8a2a1072b4a7fff, 1),
///   (0x8a2a1072b58ffff, 1),
///   (0x8a2a1072b587fff, 1),
/// ];
/// let r = libh3::k_ring_distances(0x8a2a1072b59ffff, 1)?;
/// assert_eq!(r, expected_kring_distances);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn k_ring_distances(
    origin: H3Index,
    radius: HexDistance,
) -> Result<Vec<(H3Index, HexDistance)>, H3Err> {
    let mut max = 0i64;
    match unsafe { libh3_sys::maxGridDiskSize(radius, &mut max) }.into() {
        H3Err::Success => {}
        error => return Err(error),
    }

    let mut indexes: Vec<H3Index> = vec![0; max as usize];
    let mut distances: Vec<HexDistance> = vec![0; max as usize];
    match unsafe {
        libh3_sys::gridDiskDistances(origin, radius, indexes.as_mut_ptr(), distances.as_mut_ptr())
    }
    .into()
    {
        H3Err::Success => {}
        error => return Err(error),
    }

    let mut collection =
        std::iter::zip(indexes, distances).collect::<Vec<(H3Index, HexDistance)>>();
    collection.retain(|v| v.0 != 0);
    Ok(collection)
}

/// Produce indexes within k distance of the origin index.
/// Output behavior is undefined when one of the indexes returned by this
/// function is a pentagon or is in the pentagon distortion area, H3Err::Pentagon
/// will be returned then
///
/// # Arguments
///
/// * `origin` - The center of the ring.
/// * `k_ring` - 0 is defined as the origin index, k-ring 1 is defined as k-ring 0 and
///   all neighboring indexes, and so on.
///
/// Output is placed in the provided array in order of increasing distance from
/// the origin.
pub fn hex_range(origin: H3Index, k_ring: HexDistance) -> Result<Vec<H3Index>, H3Err> {
    let mut max = 0i64;
    match unsafe { libh3_sys::maxGridDiskSize(k_ring, &mut max) }.into() {
        H3Err::Success => {}
        error => return Err(error),
    }

    let mut r: Vec<H3Index> = vec![0; max as usize];
    match unsafe { libh3_sys::gridDiskUnsafe(origin, k_ring, r.as_mut_ptr()) }.into() {
        H3Err::Success => {}
        error => return Err(error),
    }
    r.retain(|v| *v != 0);
    Ok(r)
}

/// Does the same as hex_range(), but also provides the distance as second value of the
/// returned vector
pub fn hex_range_distances(
    origin: H3Index,
    k_ring: HexDistance,
) -> Result<Vec<(H3Index, HexDistance)>, H3Err> {
    let mut max = 0i64;
    match unsafe { libh3_sys::maxGridDiskSize(k_ring, &mut max) }.into() {
        H3Err::Success => {}
        error => return Err(error),
    }

    let mut indexes: Vec<H3Index> = vec![0; max as usize];
    let mut distances: Vec<HexDistance> = vec![0; max as usize];
    match unsafe {
        libh3_sys::gridDiskDistancesUnsafe(
            origin,
            k_ring,
            indexes.as_mut_ptr(),
            distances.as_mut_ptr(),
        )
    }
    .into()
    {
        H3Err::Success => {}
        error => return Err(error),
    }
    unsafe { indexes.set_len(max as usize) };
    unsafe { distances.set_len(max as usize) };

    let result = indexes
        .into_iter()
        .zip(distances)
        .filter(|v| v.0 != 0)
        .collect::<Vec<(H3Index, HexDistance)>>();
    Ok(result)
}

/// Get the grid distance between two hex indexes. This function may fail
/// to find the distance between two indexes if they are very far apart or
/// on opposite sides of a pentagon.
/// # Arguments
///
/// * `origin` - The starting H3 index
/// * `end` - The ending H3 index
///
/// ```
/// use libh3::h3_distance;
/// assert_eq!(h3_distance(0x8a2a1072b4a7fff, 0x8a2a1072b58ffff)?, 1);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn h3_distance(origin: H3Index, end: H3Index) -> Result<i64, H3Err> {
    let mut distance = 0i64;
    match unsafe { libh3_sys::gridDistance(origin, end, &mut distance) }.into() {
        H3Err::Success => Ok(distance),
        error => Err(error),
    }
}

/// Get all hexagons with centers contained in a given polygon. The polygon
/// is specified with GeoJson semantics as an array of loops. The first loop
/// is the perimeter of the polygon, and subsequent loops are
/// expected to be holes.
///
/// # Arguments
///
/// * `polygon` - The vector of polygons.
/// * `resolution` - The resolution of the generated hexagons
///
/// ```
/// use libh3::{polyfill, GeoCoord, degs_to_rads};
/// /// Some vertexes around San Francisco
/// let sf_verts = vec![
///     (0.659966917655, -2.1364398519396),
///     (0.6595011102219, -2.1359434279405),
///     (0.6583348114025, -2.1354884206045),
///     (0.6581220034068, -2.1382437718946),
///     (0.6594479998527, -2.1384597563896),
///     (0.6599990002976, -2.1376771158464),
/// ]
/// .iter()
/// .map(|v| GeoCoord::new(v.0, v.1))
/// .collect();
///
/// let h = polyfill(&mut [sf_verts], 9)?;
/// assert_eq!(h.len(), 1253);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn polyfill(
    polygon: &mut [Vec<GeoCoord>],
    resolution: Resolution,
) -> Result<Vec<H3Index>, H3Err> {
    let mut real_polygon = polygon
        .iter()
        .map(|p| p.iter().map(libh3_sys::LatLng::from).collect())
        .collect::<Vec<Vec<libh3_sys::LatLng>>>();

    let fence = libh3_sys::GeoLoop {
        numVerts: real_polygon[0].len() as i32,
        verts: real_polygon[0].as_mut_ptr(),
    };

    let mut holes = real_polygon
        .iter()
        .skip(1)
        .map(|p| libh3_sys::GeoLoop {
            numVerts: p.len() as i32,
            verts: &mut p[0].to_owned(),
        })
        .collect::<Vec<libh3_sys::GeoLoop>>();

    let p = libh3_sys::GeoPolygon {
        geoloop: fence,
        numHoles: (real_polygon.len() - 1) as i32,
        holes: holes.as_mut_ptr(),
    };

    let mut max = 0i64;
    match unsafe {
        libh3_sys::maxPolygonToCellsSizeExperimental(
            &p,
            resolution as ::std::os::raw::c_int,
            0,
            &mut max,
        )
    }
    .into()
    {
        H3Err::Success => {}
        error => return Err(error),
    }

    let mut r: Vec<H3Index> = vec![0; max as usize];
    match unsafe {
        libh3_sys::polygonToCellsExperimental(
            &p,
            resolution as ::std::os::raw::c_int,
            0,
            max,
            r.as_mut_ptr(),
        )
    }
    .into()
    {
        H3Err::Success => {}
        error => return Err(error),
    }
    r.retain(|&v| v != 0);
    Ok(r)
}

/// Returns the size of the array needed by h3ToChildren for these inputs.
///
/// # Arguments
///
/// * `h` - The index of the parent resolution.
/// * `resolution` - The resolution of the desired level.
///
/// ```
/// use libh3::max_h3_to_children_size;
/// assert_eq!(max_h3_to_children_size(0x852a1073fffffff, 6)?, 7);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn max_h3_to_children_size(h: H3Index, resolution: Resolution) -> Result<i64, H3Err> {
    let mut out = 0i64;
    match unsafe { libh3_sys::cellToChildrenSize(h, resolution as i32, &mut out) }.into() {
        H3Err::Success => Ok(out),
        error => Err(error),
    }
}

/// Returns the parent (coarser) index containing h3.
///
/// # Arguments
///
/// * `h` - The index of the child resolution.
/// * `resolution` - The resolution of the desired level.
///
/// ```
/// use libh3::h3_to_parent;
/// assert_eq!(h3_to_parent(0x8a2a1072b4a7fff, 5)?, 0x852a1073fffffff);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn h3_to_parent(h: H3Index, resolution: Resolution) -> Result<H3Index, H3Err> {
    let mut out: H3Index = 0;
    match unsafe { libh3_sys::cellToParent(h, resolution as i32, &mut out) }.into() {
        H3Err::Success => Ok(out),
        error => Err(error),
    }
}

/// Returns children indexes contained by the given index at the given resolution.
///
/// # Arguments
///
/// * `h` - The index of the child resolution.
/// * `resolution` - The resolution of the desired level.
///
/// ```
/// use libh3::h3_to_children;
/// assert_eq!(
///     h3_to_children(0x852a1073fffffff, 6)?,
///     vec![
///         0x862a10707ffffff,
///         0x862a1070fffffff,
///         0x862a10717ffffff,
///         0x862a1071fffffff,
///         0x862a10727ffffff,
///         0x862a1072fffffff,
///         0x862a10737ffffff
///     ]
/// );
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn h3_to_children(h: H3Index, resolution: Resolution) -> Result<Vec<H3Index>, H3Err> {
    let max = max_h3_to_children_size(h, resolution)?;
    let mut result: Vec<H3Index> = vec![0; max as usize];
    match unsafe { libh3_sys::cellToChildren(h, resolution as i32, result.as_mut_ptr()) }.into() {
        H3Err::Success => {}
        error => return Err(error),
    }
    Ok(result)
}

/// Number of resolution 0 H3 indexes.
///
/// ```
/// use libh3::res_0_index_count;
/// assert_eq!(res_0_index_count(), 122);
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn res_0_index_count() -> i32 {
    unsafe { libh3_sys::res0CellCount() }
}

/// All the resolution 0 H3 indexes.
///
/// ```
/// use libh3::get_res_0_indexes;
/// assert_eq!(
///     get_res_0_indexes()?.len(),
///     122
/// );
/// # Ok::<(), libh3::H3Err>(())
/// ```
pub fn get_res_0_indexes() -> Result<Vec<H3Index>, H3Err> {
    let max = res_0_index_count() as usize;
    let mut result: Vec<H3Index> = vec![0; max as usize];
    match unsafe { libh3_sys::getRes0Cells(result.as_mut_ptr()) }.into() {
        H3Err::Success => {}
        error => return Err(error),
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_libh3() -> Result<(), H3Err> {
        // Some vertexes around San Francisco
        let sf_verts = vec![
            (0.659966917655, -2.1364398519396),
            (0.6595011102219, -2.1359434279405),
            (0.6583348114025, -2.1354884206045),
            (0.6581220034068, -2.1382437718946),
            (0.6594479998527, -2.1384597563896),
            (0.6599990002976, -2.1376771158464),
        ]
        .iter()
        .map(|v| GeoCoord::new(v.0, v.1))
        .collect();

        let h = polyfill(&mut [sf_verts], 9)?;
        assert_eq!(h.len(), 1253);

        // Fill a polygon around Wellington, NZ
        // Coordinates are in GeoJSON lng, lat order.
        let wellington_verts: Vec<GeoCoord> = vec![
            (174.800937866947, -41.22501356278325),
            (174.8079721211159, -41.226732341115365),
            (174.82262997231396, -41.231639803277986),
            (174.83561105648377, -41.23873115217201),
            (174.84634815587896, -41.24769587535717),
            (174.85069634833735, -41.252194466801384),
            (174.8587207192276, -41.26264015566857),
            (174.8636809159909, -41.27410982273725),
            (174.86536017144866, -41.28610182569948),
            (174.8653611411562, -41.29165993179072),
            (174.8636858034186, -41.30364998147521),
            (174.85872883206798, -41.31511423541933),
            (174.8507068877321, -41.32555201329627),
            (174.846359294586, -41.33004662498431),
            (174.8356222320399, -41.33900220579377),
            (174.82263983163466, -41.346084796073406),
            (174.807979604528, -41.35098543989819),
            (174.80094378133927, -41.35270175860709),
            (174.78524618901284, -41.35520670405109),
            (174.76919781098724, -41.35520670405109),
            (174.75350021866086, -41.35270175860712),
            (174.7464643954721, -41.35098543989822),
            (174.7318041683653, -41.346084796073406),
            (174.71882176795995, -41.33900220579369),
            (174.7080847054138, -41.330046624984135),
            (174.70373711226773, -41.32555201329609),
            (174.69571516793187, -41.31511423541913),
            (174.69075819658127, -41.303649981474955),
            (174.68908285884382, -41.29165993179046),
            (174.68908382855136, -41.2861018256992),
            (174.69076308400918, -41.274109822737074),
            (174.69572328077246, -41.26264015566849),
            (174.70374765166264, -41.252194466801384),
            (174.70809584412103, -41.24769587535717),
            (174.71883294351622, -41.23873115217201),
            (174.731814027686, -41.231639803278014),
            (174.746471878884, -41.22673234111539),
            (174.75350613305287, -41.22501356278328),
            (174.7691998725514, -41.222504896122565),
            (174.78524412744844, -41.222504896122565),
            (174.800937866947, -41.22501356278325),
        ]
        .iter()
        .map(|v| GeoCoord::new(degs_to_rads(v.1), degs_to_rads(v.0)))
        .collect();

        let mut h = polyfill(&mut [wellington_verts], 6)?;
        assert_eq!(h.len(), 5);
        h.sort_unstable();
        assert_eq!(
            h,
            vec![
                606774924341673983,
                606774925281198079,
                606774925549633535,
                606774929307729919,
                606774929441947647
            ]
        );
        Ok(())
    }
}
