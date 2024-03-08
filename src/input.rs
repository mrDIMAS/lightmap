use fyrox_math::octree::Octree;
use nalgebra::{Vector2, Vector3};

pub struct WorldVertex {
    pub world_normal: Vector3<f32>,
    pub world_position: Vector3<f32>,
    pub second_tex_coord: Vector2<f32>,
}

pub struct Mesh {
    /// World-space vertices.
    pub vertices: Vec<WorldVertex>,
    pub triangles: Vec<[u32; 3]>,
    pub octree: Octree,
}

impl Mesh {
    pub fn new(vertices: Vec<WorldVertex>, triangles: Vec<[u32; 3]>) -> Option<Self> {
        let mut world_space_triangles = Vec::with_capacity(triangles.len() * 3);

        for triangle in triangles.iter() {
            let a = vertices.get(triangle[0] as usize)?;
            let b = vertices.get(triangle[1] as usize)?;
            let c = vertices.get(triangle[2] as usize)?;

            world_space_triangles.push([a.world_position, b.world_position, c.world_position]);
        }

        let octree = Octree::new(&world_space_triangles, 64);

        Some(Self {
            vertices,
            triangles,
            octree,
        })
    }
}
