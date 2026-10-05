# cn_pipeline_jen_v2

Copy of Jen's updated code from the Drive (jen/jen_v2, 26 Sep): main_exec.m,
aDMDc.m, make_scenario_connectivity.m, dmf_get_params.m, balloonWindkessel.m.
The other files in jen_v2 are identical to GitHub (before Elena's 28 Sep
change) and were taken from there. Nothing here is changed.

nur_analysis/nur_jen_v2_check.m – diagnostic of this code:
  1. dt = 1 vs dt = 0.01
  2. the SVD step inside aDMDc
  3. recovery on random networks: aDMDc as is vs fixed SVD (units / areas)

Run: cd nur_analysis; nur_jen_v2_check
