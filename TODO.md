# TODO's
---
## Feature requests
   - z-directed (conventional) gradients [B]
   - custom target fields 
      - sawtooth
   - Constraints (length / radius / linearity)
   - Optimisation
      - exhaustive search fix one optimise the other
      - Length
      - Power
      - Inductance
   - Power dissipation computation [T]
   - Planar gradients
   - Load custom target field (z / Gamma(z))
   - Wire export (format?)
   - Field export (format?)
   - Minimum inductance gradients
   - Warn user when inputs lead to instable / non buildable gradients
   - Allow custom gridding (z/k/m/phi)
   - remove asserts -> move to GUI or input validation from wrapper scripts 
   - Implement non uniform DFT for grid flexibility [P]
   - Add transferfunction methods to modularise current density computation
   - Dynamically determine apodisation factor requirements
   - Determine 3D plot sizes (and gui size) based on screen size?
   - Gradient class definition with settings / generated output?

## Bug fixes
   - Avoid open wire paths
   - Robust inductance computation (and accurate!)
   - Correct for DC frequency robustly
   - Contour memory leak [T]
   - check asymptotic behaviour of Bessels in P and Q and L computation and perhaps add internal corrections when necessary?
   - parse (check) valid linearityterms
      - transverse: linearityterm even & > 3
      - longitudinal: " even & > 1
      - mexicanhat: " even & > 0 
   - Figure out a check for b < a (i.e. what radius should B be at in general?)
   - Error computation for DSV is wrong for y/(zorx) gradients (direction wrong with convention new methods?)
   - UI error, text of DSV error description doesn't update for change in DSV value.

## Code Health
   - add description to each function / module
   - Complete readme.md
   - Fix notation / coordinate system thoroughly
   - Write Manual (readthedocs.io?)
   - Speedup P and Q computation using dedicated i1e k1e k0/k1 functions
   - Remove duplicate computations from GUI/gradComp toolbox
   - how to determine DSV nr of points for accurate comparison?
