use nalgebra::{Vector2, Vector3};

#[inline]
pub fn is_point_inside_triangle(p: &Vector3<f32>, vertices: &[Vector3<f32>; 3]) -> bool {
    let ba = vertices[1] - vertices[0];
    let ca = vertices[2] - vertices[0];
    let vp = *p - vertices[0];

    let ba_dot_ba = ba.dot(&ba);
    let ca_dot_ba = ca.dot(&ba);
    let ca_dot_ca = ca.dot(&ca);

    let dot02 = ca.dot(&vp);
    let dot12 = ba.dot(&vp);

    let inv_denom = 1.0 / (ca_dot_ca * ba_dot_ba - ca_dot_ba.powi(2));

    // Calculate barycentric coordinates
    let u = (ba_dot_ba * dot02 - ca_dot_ba * dot12) * inv_denom;
    let v = (ca_dot_ca * dot12 - ca_dot_ba * dot02) * inv_denom;

    (u >= 0.0) && (v >= 0.0) && (u + v < 1.0)
}

#[inline]
pub fn barycentric_is_inside(bary: (f32, f32, f32)) -> bool {
    (bary.0 >= 0.0) && (bary.1 >= 0.0) && (bary.0 + bary.1 < 1.0)
}

#[inline]
pub fn get_barycentric_coords_2d(
    p: Vector2<f32>,
    a: Vector2<f32>,
    b: Vector2<f32>,
    c: Vector2<f32>,
) -> (f32, f32, f32) {
    let v0 = b - a;
    let v1 = c - a;
    let v2 = p - a;

    let d00 = v0.dot(&v0);
    let d01 = v0.dot(&v1);
    let d11 = v1.dot(&v1);
    let d20 = v2.dot(&v0);
    let d21 = v2.dot(&v1);
    let inv_denom = 1.0 / (d00 * d11 - d01.powi(2));

    let v = (d11 * d20 - d01 * d21) * inv_denom;
    let w = (d00 * d21 - d01 * d20) * inv_denom;
    let u = 1.0 - v - w;

    (u, v, w)
}

#[inline]
pub fn solve_quadratic(a: f32, b: f32, c: f32) -> Option<[f32; 2]> {
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        // No real roots
        None
    } else {
        // Dont care if quadratic equation has only one root (discriminant == 0), this is edge-case
        // which requires additional branching instructions which is not good for branch-predictor in CPU.
        let _2a = 2.0 * a;
        let discr_root = discriminant.sqrt();
        let r1 = (-b + discr_root) / _2a;
        let r2 = (-b - discr_root) / _2a;
        Some([r1, r2])
    }
}

#[inline]
pub fn barycentric_to_world(
    bary: (f32, f32, f32),
    pa: Vector3<f32>,
    pb: Vector3<f32>,
    pc: Vector3<f32>,
) -> Vector3<f32> {
    pa.scale(bary.0) + pb.scale(bary.1) + pc.scale(bary.2)
}

#[inline]
pub fn triangle_area(a: Vector3<f32>, b: Vector3<f32>, c: Vector3<f32>) -> f32 {
    (b - a).cross(&(c - a)).norm() * 0.5
}
