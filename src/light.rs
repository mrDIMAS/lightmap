use crate::EntityId;
use nalgebra::Vector3;

/// Directional light is a light source with parallel rays. Example: Sun.
pub struct DirectionalLightDefinition<Id: EntityId> {
    /// A handle of light in the scene.
    pub handle: Id,
    /// Intensity is how bright light is. Default is 1.0.
    pub intensity: f32,
    /// Direction of light rays.
    pub direction: Vector3<f32>,
    /// Color of light.
    pub color: Vector3<f32>,
}

/// Spot light is a cone light source. Example: flashlight.
pub struct SpotLightDefinition<Id: EntityId> {
    /// A handle of light in the scene.
    pub handle: Id,
    /// Intensity is how bright light is. Default is 1.0.
    pub intensity: f32,
    /// Color of light.
    pub color: Vector3<f32>,
    /// Direction vector of light.
    pub direction: Vector3<f32>,
    /// Position of light in world coordinates.
    pub position: Vector3<f32>,
    /// Distance at which light intensity decays to zero.
    pub distance: f32,
    /// Square of distance.
    pub sqr_distance: f32,
    /// Smoothstep left bound. It is ((hotspot_cone_angle + falloff_angle_delta) * 0.5).cos()
    pub edge0: f32,
    /// Smoothstep right bound. It is (hotspot_cone_angle * 0.5).cos()
    pub edge1: f32,
}

/// Point light is a spherical light source. Example: light bulb.
pub struct PointLightDefinition<Id: EntityId> {
    /// A handle of light in the scene.
    pub handle: Id,
    /// Intensity is how bright light is. Default is 1.0.
    pub intensity: f32,
    /// Position of light in world coordinates.
    pub position: Vector3<f32>,
    /// Color of light.
    pub color: Vector3<f32>,
    /// Radius of sphere at which light intensity decays to zero.
    pub radius: f32,
    /// Square of radius.
    pub sqr_radius: f32,
}

/// Light definition for lightmap rendering.
pub enum LightDefinition<Id: EntityId> {
    /// See docs of [DirectionalLightDefinition](struct.PointLightDefinition.html)
    Directional(DirectionalLightDefinition<Id>),
    /// See docs of [SpotLightDefinition](struct.SpotLightDefinition.html)
    Spot(SpotLightDefinition<Id>),
    /// See docs of [PointLightDefinition](struct.PointLightDefinition.html)
    Point(PointLightDefinition<Id>),
}

impl<Id: EntityId> LightDefinition<Id> {
    pub fn handle(&self) -> Id {
        match self {
            LightDefinition::Directional(v) => v.handle,
            LightDefinition::Spot(v) => v.handle,
            LightDefinition::Point(v) => v.handle,
        }
    }
}
