//! Module to generate lightmaps for surfaces.
//!
//! # Performance
//!
//! This is CPU lightmapper, its performance is linear with core count of your CPU.

#![forbid(unsafe_code)]

pub use fyrox_math as math;

pub mod error;
pub mod input;
pub mod light;

pub mod progress;

use crate::{
    input::{
        LightmapInputData,
        Mesh
    },
    error::LightmapGenerationError,
    light::LightDefinition,
    progress::{CancellationToken, ProgressIndicator, ProgressStage}
};
use arrayvec::ArrayVec;
use fxhash::FxHashMap;
use math::{
    octree::OctreeNode,
    ray::Ray,
    Rect,
};
use nalgebra::{Vector2, Vector3, Vector4};
use rayon::prelude::*;
use std::{ hash::Hash, };

pub trait EntityId: Hash + Copy + Clone + Eq + PartialEq + Send + Sync {}

pub trait Texture: Sized {
    fn from_rgb8(bytes: Vec<u8>) -> Self;
}

#[derive(Default, Clone, Debug)]
pub struct LightmapEntry<Id: EntityId, Tex: Texture> {
    /// Lightmap texture.
    ///
    /// TODO: Is single texture enough? There may be surfaces with huge amount of faces
    ///  which may not fit into texture, because there is hardware limit on most GPUs
    ///  up to 8192x8192 pixels.
    pub texture: Option<Tex>,
    /// List of lights that were used to generate this lightmap. This list is used for
    /// masking when applying dynamic lights for surfaces with light, it prevents double
    /// lighting.
    pub lights: Vec<Id>,
}

/// Lightmap is a texture with precomputed lighting.
#[derive(Default, Clone, Debug)]
pub struct Lightmap<Id: EntityId, Tex: Texture> {
    /// Node handle to lightmap mapping. It is used to quickly get information about
    /// lightmaps for any node in scene.
    pub map: FxHashMap<Id, Vec<LightmapEntry<Id, Tex>>>,
}

impl<Id: EntityId, Tex: Texture> Lightmap<Id, Tex> {
    pub fn new(
        data: LightmapInputData<Id>,
        texels_per_unit: u32,
        cancellation_token: CancellationToken,
        progress_indicator: ProgressIndicator,
    ) -> Result<Self, LightmapGenerationError> {
        let LightmapInputData {
            meshes,
            lights,
        } = data;

        progress_indicator.set_stage(ProgressStage::CalculatingLight, meshes.len() as u32);

        let mut map: FxHashMap<Id, Vec<LightmapEntry<Id, Tex>>> = FxHashMap::default();
        for instance in meshes.iter() {
            if cancellation_token.is_cancelled() {
                return Err(LightmapGenerationError::Cancelled);
            }

            let lightmap =
                generate_lightmap::<Id, Tex>(instance, &meshes, &lights, texels_per_unit);
            map.entry(instance.owner).or_default().push(LightmapEntry {
                texture: Some(lightmap),
                lights: lights.iter().map(|light| light.handle()).collect(),
            });

            progress_indicator.advance_progress();
        }

        Ok(Self { map })
    }
}

/// Computes total area of triangles in surface data and returns size of square
/// in which triangles can fit.
fn estimate_size<Id: EntityId>(data: &Mesh<Id>, texels_per_unit: u32) -> u32 {
    let mut area = 0.0;
    for triangle in data.triangles.iter() {
        let a = data.vertices[triangle[0] as usize].world_position;
        let b = data.vertices[triangle[1] as usize].world_position;
        let c = data.vertices[triangle[2] as usize].world_position;
        area += math::triangle_area(a, b, c);
    }
    area.sqrt().ceil() as u32 * texels_per_unit
}

/// Calculates distance attenuation for a point using given distance to the point and
/// radius of a light.
fn distance_attenuation(distance: f32, sqr_radius: f32) -> f32 {
    let attenuation = (1.0 - distance * distance / sqr_radius).clamp(0.0, 1.0);
    attenuation * attenuation
}

/// Calculates properties of pixel (world position, normal) at given position.
fn pick<Id: EntityId>(
    uv: Vector2<f32>,
    grid: &Grid,
    data: &Mesh<Id>,
    scale: f32,
) -> Option<(Vector3<f32>, Vector3<f32>)> {
    if let Some(cell) = grid.pick(uv) {
        for triangle in cell.triangles.iter().map(|&ti| &data.triangles[ti]) {
            let ia = triangle[0] as usize;
            let ib = triangle[1] as usize;
            let ic = triangle[2] as usize;

            let uv_a = data.vertices[ia].second_tex_coord;
            let uv_b = data.vertices[ib].second_tex_coord;
            let uv_c = data.vertices[ic].second_tex_coord;

            let center = (uv_a + uv_b + uv_c).scale(1.0 / 3.0);
            let to_center = (center - uv)
                .try_normalize(std::f32::EPSILON)
                .unwrap_or_default()
                .scale(scale * 0.3333333);

            let mut current_uv = uv;
            for _ in 0..3 {
                let barycentric = math::get_barycentric_coords_2d(current_uv, uv_a, uv_b, uv_c);

                if math::barycentric_is_inside(barycentric) {
                    let a = data.vertices[ia].world_position;
                    let b = data.vertices[ib].world_position;
                    let c = data.vertices[ic].world_position;

                    let na = data.vertices[ia].world_normal;
                    let nb = data.vertices[ib].world_normal;
                    let nc = data.vertices[ic].world_normal;

                    return Some((
                        math::barycentric_to_world(barycentric, a, b, c),
                        math::barycentric_to_world(barycentric, na, nb, nc),
                    ));
                }

                // Offset uv to center for conservative rasterization.
                current_uv += to_center;
            }
        }
    }
    None
}

struct GridCell {
    // List of triangle indices.
    triangles: Vec<usize>,
}

struct Grid {
    cells: Vec<GridCell>,
    size: usize,
    fsize: f32,
}

impl Grid {
    /// Creates uniform grid where each cell contains list of triangles
    /// whose second texture coordinates intersects with it.
    fn new<Id: EntityId>(data: &Mesh<Id>, size: usize) -> Self {
        let mut cells = Vec::with_capacity(size);
        let fsize = size as f32;
        for y in 0..size {
            for x in 0..size {
                let bounds =
                    Rect::new(x as f32 / fsize, y as f32 / fsize, 1.0 / fsize, 1.0 / fsize);

                let mut triangles = Vec::new();

                for (triangle_index, triangle) in data.triangles.iter().enumerate() {
                    let uv_a = data.vertices[triangle[0] as usize].second_tex_coord;
                    let uv_b = data.vertices[triangle[1] as usize].second_tex_coord;
                    let uv_c = data.vertices[triangle[2] as usize].second_tex_coord;
                    let uv_min = uv_a.inf(&uv_b).inf(&uv_c);
                    let uv_max = uv_a.sup(&uv_b).sup(&uv_c);
                    let triangle_bounds =
                        Rect::new(uv_min.x, uv_min.y, uv_max.x - uv_min.x, uv_max.y - uv_min.y);
                    if triangle_bounds.intersects(bounds) {
                        triangles.push(triangle_index);
                    }
                }

                cells.push(GridCell { triangles })
            }
        }

        Self {
            cells,
            size,
            fsize: size as f32,
        }
    }

    fn pick(&self, v: Vector2<f32>) -> Option<&GridCell> {
        let ix = (v.x * self.fsize) as usize;
        let iy = (v.y * self.fsize) as usize;
        self.cells.get(iy * self.size + ix)
    }
}

/// https://en.wikipedia.org/wiki/Lambert%27s_cosine_law
fn lambertian(light_vec: Vector3<f32>, normal: Vector3<f32>) -> f32 {
    normal.dot(&light_vec).max(0.0)
}

/// https://en.wikipedia.org/wiki/Smoothstep
fn smoothstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    let k = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    k * k * (3.0 - 2.0 * k)
}

/// Generates lightmap for given surface data with specified transform.
///
/// # Performance
///
/// This method is has linear complexity - the more complex mesh you pass, the more
/// time it will take. Required time increases drastically if you enable shadows and
/// global illumination (TODO), because in this case your data will be raytraced.
fn generate_lightmap<Id: EntityId, Tex: Texture>(
    instance: &Mesh<Id>,
    other_instances: &[Mesh<Id>],
    lights: &[LightDefinition<Id>],
    texels_per_unit: u32,
) -> Tex {
    // We have to re-generate new set of world-space vertices because UV generator
    // may add new vertices on seams.
    let atlas_size = estimate_size(instance, texels_per_unit);
    let scale = 1.0 / atlas_size as f32;
    let grid = Grid::new(instance, (atlas_size / 32).max(4) as usize);

    let mut pixels: Vec<Vector4<u8>> =
        vec![Vector4::new(0, 0, 0, 0); (atlas_size * atlas_size) as usize];

    let half_pixel = scale * 0.5;
    pixels
        .par_iter_mut()
        .enumerate()
        .for_each(|(i, pixel): (usize, &mut Vector4<u8>)| {
            let x = i as u32 % atlas_size;
            let y = i as u32 / atlas_size;

            let uv = Vector2::new(x as f32 * scale + half_pixel, y as f32 * scale + half_pixel);

            if let Some((world_position, world_normal)) = pick(uv, &grid, instance, scale) {
                let mut pixel_color = Vector3::default();
                for light in lights {
                    let (light_color, mut attenuation, light_position) = match light {
                        LightDefinition::Directional(directional) => {
                            let attenuation = directional.intensity
                                * lambertian(directional.direction, world_normal);
                            (directional.color, attenuation, Vector3::default())
                        }
                        LightDefinition::Spot(spot) => {
                            let d = spot.position - world_position;
                            let distance = d.norm();
                            let light_vec = d.scale(1.0 / distance);
                            let spot_angle_cos = light_vec.dot(&spot.direction);
                            let cone_factor = smoothstep(spot.edge0, spot.edge1, spot_angle_cos);
                            let attenuation = cone_factor
                                * spot.intensity
                                * lambertian(light_vec, world_normal)
                                * distance_attenuation(distance, spot.sqr_distance);
                            (spot.color, attenuation, spot.position)
                        }
                        LightDefinition::Point(point) => {
                            let d = point.position - world_position;
                            let distance = d.norm();
                            let light_vec = d.scale(1.0 / distance);
                            let attenuation = point.intensity
                                * lambertian(light_vec, world_normal)
                                * distance_attenuation(distance, point.sqr_radius);
                            (point.color, attenuation, point.position)
                        }
                    };
                    // Shadows
                    if attenuation >= 0.01 {
                        let mut query_buffer = ArrayVec::<usize, 64>::new();
                        let shadow_bias = 0.01;
                        let ray = Ray::from_two_points(light_position, world_position);
                        'outer_loop: for other_instance in other_instances {
                            other_instance
                                .octree
                                .ray_query_static(&ray, &mut query_buffer);
                            for &node in query_buffer.iter() {
                                match other_instance.octree.node(node) {
                                    OctreeNode::Leaf { indices, .. } => {
                                        let other_data = other_instance;
                                        for &triangle_index in indices {
                                            let triangle =
                                                &other_data.triangles[triangle_index as usize];
                                            let va = other_data.vertices[triangle[0] as usize]
                                                .world_position;
                                            let vb = other_data.vertices[triangle[1] as usize]
                                                .world_position;
                                            let vc = other_data.vertices[triangle[2] as usize]
                                                .world_position;
                                            if let Some(pt) =
                                                ray.triangle_intersection_point(&[va, vb, vc])
                                            {
                                                if ray.origin.metric_distance(&pt) + shadow_bias
                                                    < ray.dir.norm()
                                                {
                                                    attenuation = 0.0;
                                                    break 'outer_loop;
                                                }
                                            }
                                        }
                                    }
                                    OctreeNode::Branch { .. } => unreachable!(),
                                }
                            }
                        }
                    }
                    pixel_color += light_color.scale(attenuation);
                }

                *pixel = Vector4::new(
                    (pixel_color.x.clamp(0.0, 1.0) * 255.0) as u8,
                    (pixel_color.y.clamp(0.0, 1.0) * 255.0) as u8,
                    (pixel_color.z.clamp(0.0, 1.0) * 255.0) as u8,
                    255, // Indicates that this pixel was "filled"
                );
            }
        });

    // Prepare light map for bilinear filtration. This step is mandatory to prevent bleeding.
    let mut rgb_pixels: Vec<Vector3<u8>> = Vec::with_capacity((atlas_size * atlas_size) as usize);
    for y in 0..(atlas_size as i32) {
        for x in 0..(atlas_size as i32) {
            let fetch = |dx: i32, dy: i32| -> Option<Vector3<u8>> {
                pixels
                    .get(((y + dy) * (atlas_size as i32) + x + dx) as usize)
                    .and_then(|p| {
                        if p.w != 0 {
                            Some(Vector3::new(p.x, p.y, p.z))
                        } else {
                            None
                        }
                    })
            };

            let src_pixel = pixels[(y * (atlas_size as i32) + x) as usize];
            if src_pixel.w == 0 {
                // Check neighbour pixels marked as "filled" and use it as value.
                if let Some(west) = fetch(-1, 0) {
                    rgb_pixels.push(west);
                } else if let Some(east) = fetch(1, 0) {
                    rgb_pixels.push(east);
                } else if let Some(north) = fetch(0, -1) {
                    rgb_pixels.push(north);
                } else if let Some(south) = fetch(0, 1) {
                    rgb_pixels.push(south);
                } else if let Some(north_west) = fetch(-1, -1) {
                    rgb_pixels.push(north_west);
                } else if let Some(north_east) = fetch(1, -1) {
                    rgb_pixels.push(north_east);
                } else if let Some(south_east) = fetch(1, 1) {
                    rgb_pixels.push(south_east);
                } else if let Some(south_west) = fetch(-1, 1) {
                    rgb_pixels.push(south_west);
                } else {
                    rgb_pixels.push(Vector3::new(0, 0, 0));
                }
            } else {
                rgb_pixels.push(Vector3::new(src_pixel.x, src_pixel.y, src_pixel.z))
            }
        }
    }

    // Blur lightmap using simplest box filter.
    let mut bytes = Vec::with_capacity((atlas_size * atlas_size * 3) as usize);
    for y in 0..(atlas_size as i32) {
        for x in 0..(atlas_size as i32) {
            if x < 1 || y < 1 || x + 1 == atlas_size as i32 || y + 1 == atlas_size as i32 {
                bytes.extend_from_slice(
                    rgb_pixels[(y * (atlas_size as i32) + x) as usize].as_slice(),
                );
            } else {
                let fetch = |dx: i32, dy: i32| -> Vector3<i16> {
                    let u8_pixel = rgb_pixels[((y + dy) * (atlas_size as i32) + x + dx) as usize];
                    Vector3::new(u8_pixel.x as i16, u8_pixel.y as i16, u8_pixel.z as i16)
                };

                let north_west = fetch(-1, -1);
                let north = fetch(0, -1);
                let north_east = fetch(1, -1);
                let west = fetch(-1, 0);
                let center = fetch(0, 0);
                let east = fetch(1, 0);
                let south_west = fetch(-1, 1);
                let south = fetch(0, 1);
                let south_east = fetch(-1, 1);

                let sum = north_west
                    + north
                    + north_east
                    + west
                    + center
                    + east
                    + south_west
                    + south
                    + south_east;

                bytes.push((sum.x / 9).clamp(0, 255) as u8);
                bytes.push((sum.y / 9).clamp(0, 255) as u8);
                bytes.push((sum.z / 9).clamp(0, 255) as u8);
            }
        }
    }

    Texture::from_rgb8(bytes)
}
