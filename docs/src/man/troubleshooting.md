# Troubleshooting

## Common Issues

1. **FrameSystem Setup**: Ensure your `FrameSystem` has all required frames (ICRF, body-fixed) and points (Sun, Moon, central body) before creating `FrameAwareParams`.
2. **NAIF IDs**: Use correct NAIF IDs in `FrameEphemeris` (399=Earth, 10=Sun, 301=Moon, 499=Mars, 2000433=Eros, etc.).
3. **Units**: All positions in km, velocities in km/s. Gravity models expect m internally (conversion is handled automatically).
4. **Earth-only models**: Drag, magnetic field, and plasma drag require EOP data and ITRF frame.
5. **`R_Occulting`**: SRP and thermal models require an explicit occulting body radius -- no default is assumed.
6. **`frames` kwarg**: Passing `frames` enables allocation-free compiled transforms via `compile_rotation`/`compile_translation`. Multi-hop paths in the frame graph are supported; compilation failures are caught and the model transparently falls back to runtime `rotation3`/`vector3`/`vector6` lookups.
7. **Keplerian gravity**: Requires an explicit `μ` parameter -- no default is assumed.
