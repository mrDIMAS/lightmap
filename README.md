# lightmap

Light maps generator.

## Example

```rust
use lightmap::{
    input::{Mesh, WorldVertex},
    light::{LightDefinition, PointLightDefinition},
    LightMap,
};
use nalgebra::{Vector2, Vector3};

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Vertex {
    pub position: Vector3<f32>,
    pub tex_coord: Vector2<f32>,
}

impl Vertex {
    fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            position: Vector3::new(x, y, z),
            tex_coord: Default::default(),
        }
    }
}

// Create cube geometry.
let mut vertices = vec![
    Vertex::new(-0.5, -0.5, 0.5),
    Vertex::new(-0.5, 0.5, 0.5),
    Vertex::new(0.5, 0.5, 0.5),
    Vertex::new(0.5, -0.5, 0.5),
    Vertex::new(-0.5, -0.5, -0.5),
    Vertex::new(-0.5, 0.5, -0.5),
    Vertex::new(0.5, 0.5, -0.5),
    Vertex::new(0.5, -0.5, -0.5),
];
let mut triangles = vec![
    // Front
    [2, 1, 0],
    [3, 2, 0],
    // Back
    [4, 5, 6],
    [4, 6, 7],
    // Right
    [7, 6, 2],
    [2, 3, 7],
    // Left
    [0, 1, 5],
    [0, 5, 4],
    // Top
    [5, 1, 2],
    [5, 2, 6],
    // Bottom
    [3, 0, 4],
    [7, 3, 4],
];

// Step 1. Generate second texture coordinates first.
let patch = uvgen::generate_uvs(
    vertices.iter().map(|v| v.position),
    triangles.iter().cloned(),
    0.005,
)
.unwrap();

// Step 2. Clone vertices on seams.
for &vertex_index in &patch.additional_vertices {
    let vertex = vertices[vertex_index as usize];
    vertices.push(vertex);
}

// Step 3. Assign generated second texture coordinates.
for (vertex, tex_coord) in vertices.iter_mut().zip(&patch.second_tex_coords) {
    vertex.tex_coord = *tex_coord;
}

// Step 4. Replace topology of the mesh using new one, that was produced when generating UVs.
triangles = patch.triangles;

// Step 5. Generate light map.
let lights = [LightDefinition::Point(PointLightDefinition {
    intensity: 1.0,
    color: Vector3::new(1.0, 1.0, 1.0),
    radius: 2.0,
    position: Vector3::new(0.0, 2.0, 0.0),
    sqr_radius: 4.0,
})];

let mut mesh = Mesh::new(
    vertices
        .iter()
        .map(|v| WorldVertex {
            world_normal: Default::default(),
            world_position: v.position,
            second_tex_coord: v.tex_coord,
        })
        .collect(),
    triangles,
)
.unwrap();

// Calculate normals.
for triangle in mesh.triangles.iter() {
    let pos_a = mesh.vertices[triangle[0] as usize].world_position;
    let pos_b = mesh.vertices[triangle[1] as usize].world_position;
    let pos_c = mesh.vertices[triangle[2] as usize].world_position;

    let normal = (pos_c - pos_a)
        .cross(&(pos_c - pos_b))
        .try_normalize(f32::EPSILON)
        .unwrap_or_default();

    mesh.vertices[triangle[0] as usize].world_normal = normal;
    mesh.vertices[triangle[1] as usize].world_normal = normal;
    mesh.vertices[triangle[2] as usize].world_normal = normal;
}

let meshes = vec![mesh];

let light_map = LightMap::new(&meshes[0], &meshes, &lights, 64);

dbg!(light_map);
```
