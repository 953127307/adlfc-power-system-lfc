# ADLFC for nonlinear power-system load frequency control

MATLAB simulation code accompanying **A Novel Adaptive Data-Driven Load Frequency Control for Unknown Nonlinear Power Systems**.

This snapshot reproduces the numerical experiments associated with manuscript Figs. 3-13 and Tables III-V. The simulation kernels are copied without numerical changes from the versions used for those figures. Historical project folders are not needed.

[中文说明](README_zh.md)

## Run

Tested with **MATLAB R2023b** on Windows. The simulations, controller reconstruction, metrics, and plots use MATLAB functions; Simulink and external Python packages are not required.

Clone this repository, change MATLAB's current folder to its root, and run:

```matlab
run_all
```

This runs six full-resolution cases, exports the three paper tables, checks the reproduced trajectories against saved reference samples, and generates eleven PNG figures in `results/figures/`.

Individual stages:

```matlab
run_experiment_a       % Section IV-A: three noise levels, 600 s
run_experiment_b       % Section IV-B: three control methods, 420 s
compute_paper_metrics  % Tables III-V and system-wide ACE peak reductions
verify_results        % Numerical reference and table checks
plot_results          % Figures from regenerated results
```

The default simulations use the exact frozen controller weights used for the manuscript figures. To repeat the offline data collection and Algorithm 1 synthesis separately:

```matlab
[w, certificate] = compute_adlfc_weights;
```

The regenerated weights and recorded data are saved to `results/controller_regenerated.mat`; this command does not replace the frozen inputs.

## Experiments and numerical settings

| Setting | Experiment A | Experiment B |
|---|---|---|
| Duration | 600 s | 420 s |
| Integration | Forward Euler, 0.001 s | Forward Euler, 0.001 s |
| Initial plant states | Zero, see configuration | Recorded near-equilibrium states, see configuration |
| State-equation disturbance bound | 1e-10, 0.05, 0.1 for Cases 1-3 | 1e-10 for all methods |
| Governor derivative limit | +/-2 p.u./s | +/-2 p.u./s |
| Control input amplitude clipping | None | None |
| Filter time constants | 0.01, 0.01, 0.02 s | 0.01, 0.01, 0.02 s |
| Adaptive gains r2, r3, r4 | 0.5, 0.5, 10 | 0.5, 0.5, 10 |
| Leakage gains delta2, delta3, delta4 | 0.01, 0.001, 0.01 | 0.01, 0.001, 0.01 |

The nominal noise-free cases retain the simulation's numerical bound `1e-10`; they are not rerun with a different bound. State-equation noise uses `rng(1)` in each simulation kernel.

Experiment A uses load changes at 60, 200, 360, and 500 s. Experiment B uses a common baseline `[0.020, 0.028, 0.024, 0.026]` p.u. and the following cumulative load additions:

- Area 1: +0.0125 p.u. step at 60 s.
- Area 2: multisine from 200 s, frequencies `[0.01, 0.05, 0.2]` Hz and amplitudes `[0.003, 0.0015, 0.0005]` p.u.
- Area 3: Ornstein-Uhlenbeck load from 360 s, correlation time 5 s and stationary standard deviation parameter 0.0025 p.u.
- Area 4: constant baseline load.

The multisine and random loads use a 5 s cosine ramp. The OU generator uses an `mt19937ar` stream with seed `20261027` and starts at zero. The three controllers share the same load at every integration step. `src/build_mixed_load.m` regenerates the stored load trajectory and the B entry point checks equality before simulating.

B controller order is ADLFC, ACE-SMC with theta adaptation, and integral ACE feedback. The baseline implementation is in `src/B_mfile_sim.m`: SMC uses `W=500`, `T_theta=0.1`, `M1=2`, `M2=0.1`, `M3=0.01`, `M4=-(M2+M3)`, and `M5=0.0217*M1*K_p/mu_p`. Integral ACE uses gain `0.2`. SMC auxiliary states, theta, and the control/ACE integrators start at zero. Exact plant and ADLFC filter initial states are in `config/frozen_settings.json` and its MAT equivalent.

## ACE and reported metrics

For area i, the implemented ACE is

```text
ACE_i = sum_l B(i,l) * (P_f,l - P_f,l*) + b_i * gamma_i
b_i   = 1/R_i + 1/K_p,i
```

Here `B` is the incidence matrix, line powers and ACE are in p.u., and `gamma_i` is in rad/s. The implementation uses this angular-frequency variable directly, without an additional Hz conversion. The four bias coefficients are approximately `[0.4083333333, 0.3792592593, 0.3933110368, 0.3655816757]` p.u./(rad/s). The system-wide quantity is `E_ACE(t) = max_i abs(ACE_i(t))`.

- RoCoF: maximum absolute forward-difference frequency derivative over the first 0.5 s after each event.
- Nadir: minimum frequency deviation over the first 10 s after each load-increase event.
- ACE settling time: time from the event to the first re-entry followed by at least one full second within the band. The zero-centered half-width is 1% of the largest absolute ACE across all four areas over that A-case run. Cases already inside the band throughout an event interval report 0; lack of such a re-entry before the next event/end reports Inf.
- System ACE peak comparison: inclusive 30 s windows after 60, 200, and 360 s; reduction is relative to integral ACE.

The copied metrics helper also computes historical integral-error fields internally. The paper-facing export uses only RoCoF, Nadir, and ACE settling time.

## Data and outputs

- `data/controller_design.mat`: frozen weights, offline sampled matrices, construction variables, and diagnostic certificate. Offline data collection uses a 1 s run sampled every 0.1 s, including both endpoints (11 samples). Diagnostic disturbance/model information is not an input to the weight equations.
- `config/frozen_settings.mat` and `.json`: parameters and initial states read from the figure-generating runs.
- `data/mixed_load_profile.mat`: shared full-resolution B load trajectory.
- `reference/reference_traces.mat`: compact reference samples, including regular samples, event neighbors, and signal extrema, from the original six runs.
- `reference/paper_tables.json`: the reported table values used for checks at manuscript precision.
- `docs/source_provenance.json`: copied-file origins and hashes.

Generated full-resolution results are stored under ignored `results/`. Run `verify_results` after simulation to create `results/verification.json`. Plotting uses the regenerated signals and preserves sampled extrema; the repository's simplified figure layout is not the manuscript's typesetting layout.

## Validation snapshot

The initial public snapshot was checked in MATLAB R2023b: all six regenerated full-resolution time, state, input, line-power, and ACE arrays matched the original figure-generating runs exactly (maximum absolute difference 0). All three paper tables matched at reported precision. Offline weight reconstruction also matched the frozen weights exactly. The checked snapshot is recorded in `docs/validation_summary.json`.

## Method references

The ACE-SMC comparison follows the implementation used in the accompanying manuscript, based on S. Trip, M. Cucuzzella, C. De Persis, A. Ferrara, and J. M. A. Scherpen, "Robust load frequency control of nonlinear power networks," International Journal of Control, 93(2), 346-359, 2020, DOI: 10.1080/00207179.2018.1557338. The integral ACE comparison is associated with P. Kundur, Power System Stability and Control, McGraw-Hill, 1994. The precise implemented laws and settings are included here rather than inferred from those references.

## License

No license has been added to this repository.
