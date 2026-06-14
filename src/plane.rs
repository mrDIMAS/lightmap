use nalgebra::Vector3;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Plane {
    pub normal: Vector3<f32>,
    pub d: f32,
}

impl Default for Plane {
    #[inline]
    fn default() -> Self {
        Plane {
            normal: Vector3::new(0.0, 1.0, 0.0),
            d: 0.0,
        }
    }
}

impl Plane {
    /// Creates plane from a point and normal vector at that point.
    /// May fail if normal is degenerated vector.
    #[inline]
    pub fn from_normal_and_point(normal: &Vector3<f32>, point: &Vector3<f32>) -> Option<Self> {
        normal
            .try_normalize(f32::EPSILON)
            .map(|normalized_normal| Self {
                normal: normalized_normal,
                d: -point.dot(&normalized_normal),
            })
    }

    /// Tries to create a plane from three points (triangle). May fail if the triangle is degenerated
    /// (collapsed into a point or a line).
    #[inline]
    pub fn from_triangle(a: &Vector3<f32>, b: &Vector3<f32>, c: &Vector3<f32>) -> Option<Self> {
        let normal = (b - a).cross(&(c - a));
        Self::from_normal_and_point(&normal, a)
    }

    /// Creates plane using coefficients of plane equation Ax + By + Cz + D = 0
    /// May fail if length of normal vector is zero (normal is degenerated vector).
    #[inline]
    pub fn from_abcd(a: f32, b: f32, c: f32, d: f32) -> Option<Self> {
        let normal = Vector3::new(a, b, c);
        let len = normal.norm();
        if len == 0.0 {
            None
        } else {
            let coeff = 1.0 / len;
            Some(Self {
                normal: normal.scale(coeff),
                d: d * coeff,
            })
        }
    }

    #[inline]
    pub fn dot(&self, point: &Vector3<f32>) -> f32 {
        self.normal.dot(point) + self.d
    }

    #[inline]
    pub fn distance(&self, point: &Vector3<f32>) -> f32 {
        self.dot(point).abs()
    }

    /// Projects the given point onto the plane along the normal vector of the plane.
    #[inline]
    pub fn project(&self, point: &Vector3<f32>) -> Vector3<f32> {
        point - self.normal.scale(self.normal.dot(point) + self.d)
    }

    /// <http://geomalgorithms.com/a05-_intersect-1.html>
    pub fn intersection_point(&self, b: &Plane, c: &Plane) -> Vector3<f32> {
        let f = -1.0 / self.normal.dot(&b.normal.cross(&c.normal));

        let v1 = b.normal.cross(&c.normal).scale(self.d);
        let v2 = c.normal.cross(&self.normal).scale(b.d);
        let v3 = self.normal.cross(&b.normal).scale(c.d);

        (v1 + v2 + v3).scale(f)
    }
}
