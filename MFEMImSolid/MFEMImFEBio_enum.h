#pragma once

namespace febio_im_solid {

// Type of solid scalar field being output
// Stress: derived from Cauchy stress
// Strain: derived from Green-Lagrange or appropriate strain tensor
enum class SolidFieldKind {
    Stress = 0,
    Strain = 1
};

} // namespace febio_im_solid
