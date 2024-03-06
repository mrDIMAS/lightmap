use crate::light::LightDefinition;
use crate::EntityId;
use fyrox_math::octree::Octree;
use nalgebra::{Vector2, Vector3};

pub struct WorldVertex {
    pub world_normal: Vector3<f32>,
    pub world_position: Vector3<f32>,
    pub second_tex_coord: Vector2<f32>,
}

pub struct Mesh<Id: EntityId> {
    pub owner: Id,
    /// World-space vertices.
    pub vertices: Vec<WorldVertex>,
    pub triangles: Vec<[u32; 3]>,
    pub octree: Octree,
}

/// Data set required to generate a lightmap. It could be produced from a scene using [`LightmapInputData::from_scene`] method.
/// It is used to split preparation step from the actual lightmap generation; to be able to put heavy generation in a separate
/// thread.
pub struct LightmapInputData<Id: EntityId> {
    pub meshes: Vec<Mesh<Id>>,
    pub lights: Vec<LightDefinition<Id>>,
}
