# Conceptual cold-pool propagation and source inversion

This separate study uses **only the original 197 five-minute cold-pool candidates**, with 188 Eq. (2)-supported and 9 Eq. (3)-supported episodes. The 1-minute and modified-monotonicity catalogs are not used or altered.

## Deliverables

- [ZIP containing all results and code](BAMNET_Model_Propagation_All_Results.zip).
- [English meeting PDF](Propagation_Meeting_EN.pdf), [Indonesian PDF](Propagation_Meeting_ID.pdf), 12 pages each.
- [Results workbook](Propagation_Results.xlsx).
- [Three-stage synthetic animation](animations/three_models.gif).
- [Conditional real-data animation](animations/real_conditional_example.gif).
- [Model implementation](propagation.py), [study runner](run_study.py), [validation and reports](finish_study.py), [unit tests](test_propagation.py).
- `inputs/`, `tables/`, `plots/`: the 197-event input, coordinates, associations, scenario estimates, alternative solutions, uncertainty experiments, and figures.

On GitHub, open a file and click **Download raw file**.

## Models and inverse problem

1. **Circle, no background wind:** distance from source to station equals radial speed times elapsed time. Unknowns: source x/y, release time, and radial speed c.
2. **Advected circle:** center(t)=source+U(t−t0), radius=c(t−t0). U is prescribed translation; measured surface wind is a proxy, not proven cold-pool advection.
3. **Advected ellipse:** major and minor semi-axes grow at c√r and c/√r; r≥1 is the axis ratio, and theta the orientation. r=1 reduces exactly to the circle. The arrival-time quadratic is solved for first positive passage.

Coordinates are horizontal km, times minutes, velocities m/s, and angles counterclockwise from east. These are constant-speed kinematic models, not fluid-dynamical or terrain-following simulations.

Bounded nonlinear least squares estimates a source from point arrival times, using multiple starting guesses. Multi-start is not proof of global uniqueness. A free circle with prescribed wind has four unknowns; four arrivals leave zero residual degrees of freedom. A free ellipse has six unknowns and cannot be identified from five single-pass arrivals. Unknown U adds another two parameters.

## Controlled experiments

Nine spatially distributed synthetic stations observe a source at (−12,−8) km, release time −60 minutes, and c=6 m/s. Advected cases use U=(2,1) m/s; the ellipse has r=1.6 and theta=30°. All three noiseless inverse problems recover the source to numerical precision (<0.001 km and <10⁻⁵ minute RMSE).

Thirty-five independent uniform ±2.5-minute perturbations per model give median source errors of **0.421, 0.542, and 0.537 km** for the no-wind circle, advected circle, and advected ellipse. Maximum errors are 0.823, 1.064, and 1.258 km. These are conditional synthetic timing sensitivities with correct model and prescribed wind, **not real-source confidence intervals**. Actual onset errors may exceed the optimistic half-bin perturbations.

Repeating the synthetic ellipse inversion on the actual five BAMNET coordinates produces multiple near-exact solutions with different origins and shapes. The Jacobian rank cannot exceed five for six unknowns; an excellent fit does not establish an identifiable source.

## Actual geometry and association

Station coordinates come from NetCDF metadata and are consistent across available months. A local equirectangular projection is used. The network's minor/major spatial spread ratio is **0.108**, with elevation ranging from **27 to 945 m ASL**; terrain is ignored by the model.

The 197 detections are grouped into temporal **hypotheses**, not independently verified common cold pools. A bounded greedy grouping uses at most 180 minutes from the first detection; sensitivity uses 60/90/120 minutes. At 180 minutes, there are 143 groups, 40 with at least two sites, **11 with at least three**, **2 with four**, and none with five. Repeated-station ambiguous groups are flagged and not source-fitted.

At 60/90/120 minutes, the numbers of groups with at least three sites are 9/10/11; counts with at least four sites are 0/1/1. Timing proximity and similar cooling do not prove a shared physical origin. No claim is made that 143 groups are 143 distinct physical cold pools.

## Conditional real-data source estimates

Main scenarios prescribe c=5 m/s and fit only x/y/release time; sensitivity uses 3 and 8 m/s. U is the component-wise median across station pre-event means, using −40..−20 minutes relative to each local onset, with at least four of five samples required. A constant surface-wind proxy is not independently validated advection.

The ellipse scenarios prescribe r=1.5, with orientation along U (or the network axis for weak U); a perpendicular orientation is also explored. These shapes are assumptions, not measured ellipse parameters. Real free-ellipse inversion is deliberately not claimed.

Spatial bounds are ±100 km from the local projection origin. Release time is limited to the six hours before the earliest detection. When speed is fitted, it is bounded to 1–30 m/s and constrained to subcritical advection for the shape; fixed-speed supercritical scenarios are explicitly skipped. All bounds are modeling priors. Boundary-hitting estimates and lack of residual degrees of freedom are flagged. Domain sensitivity ±50/100/150 km is supplied for the example case.

### Example G059: 28 June 2025

Observed onsets are BAMNET01 05:10, BAMNET03 05:40, BAMNET04 06:05, BAMNET05 06:20 UTC. Prescribed background U=(0.646,1.772) m/s. With c fixed at 5 m/s:

| Scenario | Source longitude | Source latitude | Release time UTC, rounded | Arrival RMSE |
|---|---:|---:|---|---:|
| No-wind circle | 100.179°E | −0.626° | 04:29 | 2.56 min |
| Advected circle | 100.102°E | −0.692° | 04:22 | 5.35 min |
| Prescribed ellipse | 100.037°E | −0.756° | 03:57 | 5.35 min |

These are **conditional hypotheses**, not confirmed geographic origins. A smaller residual for the no-wind model does not prove zero environmental wind; wind-proxy errors, association errors, and model inadequacy remain possible. Freeing c can make the residual almost zero with four data and four parameters—this is not independent validation. Another four-station group, G043 on 21 April, has larger residuals and prior-boundary solutions.

The source-loss map fixes c=5 m/s and optimizes release time at each grid point. Its white contour is best-grid RMSE+2.5 minutes, **not a 95% confidence boundary**. Alternative shapes, speeds, wind proxies, association choices, and priors contribute uncertainty beyond this map. Numerical solver digits in exported coordinates/times must not be interpreted as measurement precision. Catalog thermal onsets need not mark exactly the same front phase at every station; local responses and temperature–gust timing offsets can bias inversion. Fitted release time is an effective backward extrapolation, not a directly observed time of convective initiation.

## Animations and next steps

In the synthetic animation, green stations have been reached by the **model**. In the real conditional animation, green means an onset has been **observed in the catalog**; a model front may not align exactly. The static example also plots measured temperature traces and observed versus predicted arrivals.

Next: independently validate shared-event associations, constrain propagation speed and translation wind, add off-axis measurements, and use radar/satellite/terrain evidence. A kinematic expansion center need not equal the exact convective initiation or rainfall location.

## Reproduction and checks

Use Python 3.12 and `requirements.txt`. The Indonesian README includes full Windows commands. Run `run_study.py` with `--data-root`, `--catalog inputs/coldpool_197_5min.csv`, and `--output .`; then `finish_study.py --data-root ...`.

Seven unit tests pass. Additional checks verify the 197-event input, nine Eq. (3) episodes, exactly-once group membership, synthetic recovery, ellipse rank deficiency, real arrival predictions, GIF frame counts, and workbook rows. See `validation_report.json`. Earlier detection outputs and source NetCDF files are unchanged.
