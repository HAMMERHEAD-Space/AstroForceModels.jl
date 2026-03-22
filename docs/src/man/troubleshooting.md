# Troubleshooting

## Common Issues

1. **FrameSystem Setup**: Ensure your `FrameSystem` has all required frames (ICRF, body-fixed) and points (Sun, Moon, central body) before creating `FrameAwareParams`.
2. **NAIF IDs**: Use correct NAIF IDs in `FrameEphemeris` (399=Earth, 10=Sun, 301=Moon, 499=Mars, 2000433=Eros, etc.).
3. **Units**: All positions in km, velocities in km/s. Gravity models expect m internally (conversion is handled automatically).
4. **Earth-only models**: Drag, magnetic field, and plasma drag require EOP data and ITRF frame.
5. **`R_Occulting`**: SRP and thermal models require an explicit occulting body radius -- no default is assumed.
6. **`frames` kwarg**: When passing `frames` to enable compiled transforms, the frame pair must be a direct parent-child connection in the frame graph. All standard setups (ICRF->ITRF, ICRF->IAU\_MARS, ICRF->ErosFixed) satisfy this.
7. **Keplerian gravity**: Requires an explicit `μ` parameter -- no default is assumed.
