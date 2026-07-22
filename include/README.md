# Source Code

This directory contains the C implementations of the algorithms developed during the dissertation.

## Files

| File | Description |
|------|-------------|
| `foc.c` | Complete implementation of the Field-Oriented Control (FOC) algorithm, including current control, coordinate transformations, decoupling, and Space Vector PWM generation. |
| `Estimador_B.c` | Experimental viscous friction coefficient (B) estimation algorithm. |
| `Estimador_J.c` | Rotor inertia (J) estimation algorithm. |
| `Estimador_LqLd.c` | Estimation routine for the d- and q-axis inductances. |

## Notes

The controller implementation was originally developed and validated in PSIM before being prepared for deployment on ARM-based embedded systems.

The current version is provided as a reference implementation accompanying the MSc dissertation.
