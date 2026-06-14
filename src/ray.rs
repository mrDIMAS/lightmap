use crate::plane::Plane;
use crate::utils::{is_point_inside_triangle, solve_quadratic};
use nalgebra::Vector3;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Ray {
    pub origin: Vector3<f32>,
    pub dir: Vector3<f32>,
}

impl Default for Ray {
    #[inline]
    fn default() -> Self {
        Ray {
            origin: Vector3::new(0.0, 0.0, 0.0),
            dir: Vector3::new(0.0, 0.0, 1.0),
        }
    }
}

/// Pair of ray equation parameters.
#[derive(Clone, Debug, Copy)]
pub struct IntersectionResult {
    pub min: f32,
    pub max: f32,
}

impl IntersectionResult {
    #[inline]
    pub fn from_slice(roots: &[f32]) -> Self {
        let mut min = f32::MAX;
        let mut max = -f32::MAX;
        for n in roots {
            min = min.min(*n);
            max = max.max(*n);
        }
        Self { min, max }
    }

    #[inline]
    pub fn from_set(results: &[Option<IntersectionResult>]) -> Option<Self> {
        let mut result = None;
        for v in results {
            match result {
                None => result = *v,
                Some(ref mut result) => {
                    if let Some(v) = v {
                        result.merge(v.min);
                        result.merge(v.max);
                    }
                }
            }
        }
        result
    }

    /// Updates min and max ray equation parameters according to a new parameter -
    /// expands range if `param` was outside of that range.
    #[inline]
    pub fn merge(&mut self, param: f32) {
        if param < self.min {
            self.min = param;
        }
        if param > self.max {
            self.max = param;
        }
    }

    #[inline]
    pub fn merge_slice(&mut self, params: &[f32]) {
        for param in params {
            self.merge(*param)
        }
    }
}

impl Ray {
    /// Creates ray from two points. May fail if begin == end.
    #[inline]
    pub fn from_two_points(begin: Vector3<f32>, end: Vector3<f32>) -> Self {
        Ray {
            origin: begin,
            dir: end - begin,
        }
    }

    #[inline]
    pub fn new(origin: Vector3<f32>, dir: Vector3<f32>) -> Self {
        Self { origin, dir }
    }

    /// Checks intersection with sphere. Returns two intersection points or none
    /// if there was no intersection.
    #[inline]
    pub fn sphere_intersection_points(
        &self,
        position: &Vector3<f32>,
        radius: f32,
    ) -> Option<[Vector3<f32>; 2]> {
        self.try_eval_points(self.sphere_intersection(position, radius))
    }

    #[inline]
    pub fn sphere_intersection(
        &self,
        position: &Vector3<f32>,
        radius: f32,
    ) -> Option<IntersectionResult> {
        let d = self.origin - *position;
        let a = self.dir.dot(&self.dir);
        let b = 2.0 * self.dir.dot(&d);
        let c = d.dot(&d) - radius * radius;
        solve_quadratic(a, b, c).map(|roots| IntersectionResult::from_slice(&roots))
    }

    /// Checks intersection with sphere.
    #[inline]
    pub fn is_intersect_sphere(&self, position: &Vector3<f32>, radius: f32) -> bool {
        let d = self.origin - position;
        let a = self.dir.dot(&self.dir);
        let b = 2.0 * self.dir.dot(&d);
        let c = d.dot(&d) - radius * radius;
        let discriminant = b * b - 4.0 * a * c;
        discriminant >= 0.0
    }

    /// Returns t factor (at pt=o+d*t equation) for projection of given point at ray
    #[inline]
    pub fn project_point(&self, point: &Vector3<f32>) -> f32 {
        (point - self.origin).dot(&self.dir) / self.dir.norm_squared()
    }

    /// Returns point on ray which defined by pt=o+d*t equation.
    #[inline]
    pub fn get_point(&self, t: f32) -> Vector3<f32> {
        self.origin + self.dir.scale(t)
    }

    #[inline]
    pub fn box_intersection(
        &self,
        min: &Vector3<f32>,
        max: &Vector3<f32>,
    ) -> Option<IntersectionResult> {
        let (mut tmin, mut tmax) = if self.dir.x >= 0.0 {
            (
                (min.x - self.origin.x) / self.dir.x,
                (max.x - self.origin.x) / self.dir.x,
            )
        } else {
            (
                (max.x - self.origin.x) / self.dir.x,
                (min.x - self.origin.x) / self.dir.x,
            )
        };

        let (tymin, tymax) = if self.dir.y >= 0.0 {
            (
                (min.y - self.origin.y) / self.dir.y,
                (max.y - self.origin.y) / self.dir.y,
            )
        } else {
            (
                (max.y - self.origin.y) / self.dir.y,
                (min.y - self.origin.y) / self.dir.y,
            )
        };

        if tmin > tymax || (tymin > tmax) {
            return None;
        }
        if tymin > tmin {
            tmin = tymin;
        }
        if tymax < tmax {
            tmax = tymax;
        }
        let (tzmin, tzmax) = if self.dir.z >= 0.0 {
            (
                (min.z - self.origin.z) / self.dir.z,
                (max.z - self.origin.z) / self.dir.z,
            )
        } else {
            (
                (max.z - self.origin.z) / self.dir.z,
                (min.z - self.origin.z) / self.dir.z,
            )
        };

        if (tmin > tzmax) || (tzmin > tmax) {
            return None;
        }
        if tzmin > tmin {
            tmin = tzmin;
        }
        if tzmax < tmax {
            tmax = tzmax;
        }
        if tmin <= 1.0 && tmax >= 0.0 {
            Some(IntersectionResult {
                min: tmin,
                max: tmax,
            })
        } else {
            None
        }
    }

    /// Solves plane equation in order to find ray equation parameter.
    /// There is no intersection if result < 0.
    #[inline]
    pub fn plane_intersection(&self, plane: &Plane) -> f32 {
        let u = -(self.origin.dot(&plane.normal) + plane.d);
        let v = self.dir.dot(&plane.normal);
        u / v
    }

    #[inline]
    pub fn plane_intersection_point(&self, plane: &Plane) -> Option<Vector3<f32>> {
        let t = self.plane_intersection(plane);
        if !(0.0..=1.0).contains(&t) {
            None
        } else {
            Some(self.get_point(t))
        }
    }

    #[inline]
    pub fn triangle_intersection_point(
        &self,
        vertices: &[Vector3<f32>; 3],
    ) -> Option<Vector3<f32>> {
        let ba = vertices[1] - vertices[0];
        let ca = vertices[2] - vertices[0];
        let plane = Plane::from_normal_and_point(&ba.cross(&ca), &vertices[0])?;

        if let Some(point) = self.plane_intersection_point(&plane) {
            if is_point_inside_triangle(&point, vertices) {
                return Some(point);
            }
        }
        None
    }

    #[inline]
    pub fn try_eval_points(&self, result: Option<IntersectionResult>) -> Option<[Vector3<f32>; 2]> {
        match result {
            None => None,
            Some(result) => {
                let a = if result.min >= 0.0 && result.min <= 1.0 {
                    Some(self.get_point(result.min))
                } else {
                    None
                };

                let b = if result.max >= 0.0 && result.max <= 1.0 {
                    Some(self.get_point(result.max))
                } else {
                    None
                };

                match a {
                    None => b.map(|b| [b, b]),
                    Some(a) => match b {
                        None => Some([a, a]),
                        Some(b) => Some([a, b]),
                    },
                }
            }
        }
    }
}
