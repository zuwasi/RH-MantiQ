#!/usr/bin/env wolframscript
projectDir = DirectoryName[$InputFileName];
SetDirectory[projectDir];
exportDir = FileNameJoin[{projectDir, "exports"}];
If[! DirectoryQ[exportDir], CreateDirectory[exportDir]];

Get[FileNameJoin[{projectDir, "src", "RHVerification.wl"}]];
LaunchKernels[];

scale = Environment["RH_SCALE"];
scale = If[StringQ[scale] && scale == "demo", "demo", "full"];

config = If[scale == "demo",
  <|"lagariasN" -> 10^6, "zeroK" -> 2000, "psiX" -> 500, "psiK" -> 200, "apQ" -> 200|>,
  <|"lagariasN" -> 10^8, "zeroK" -> 20000, "psiX" -> 2000, "psiK" -> 1000, "apQ" -> 1000|>
];
Print["SCALE: ", scale, "  config: ", config];

machine = MachineProfile[];
Export[FileNameJoin[{exportDir, "machine.json"}], machine, "JSON"];
Print["Machine: ", machine["WolframVersion"], " x ", machine["ProcessorCount"], " cores"];

(* ---- 1. Lagarias RH-equivalence ---- *)
Print["\n=== Lagarias criterion up to n = ", config["lagariasN"], " ==="];
lag = LagariasCheck[config["lagariasN"]];
lagExport = KeyDrop[lag, "tightest"];
AppendTo[lagExport, "tightest" -> Normal[lag["tightest"]]];
Export[FileNameJoin[{exportDir, "lagarias_summary.json"}], lagExport, "JSON"];
Export[FileNameJoin[{exportDir, "lagarias_tight.csv"}], lag["tightest"], "CSV"];
Print["violations: ", lag["violations"], "  recheckNeeded: ", lag["recheckNeeded"],
  "  maxRatio (n>1): ", N[lag["maxRatio"], 6], "  elapsed: ", N[lag["elapsedSeconds"], 4], "s"];

figTight = Labeled[
  BarChart[N[lag["tightest"][[All, "ratio"]]],
   ChartLabels -> Style[#, 12] & /@ (lag["tightest"][[All, "n"]]),
   AxesLabel -> {"n", "sigma(n)/bound(n)"},
   PlotLabel -> "Lagarias ratio: " <> ToString[config["lagariasN"]] <> " range, tightest n (n>1)"],
  "RH equivalence: sigma(n) <= H_n + E^H_n Log[H_n]", Top];
Export[FileNameJoin[{exportDir, "fig_lagarias_tight.png"}], figTight, ImageResolution -> 200];

(* ---- 2. Zeros on the critical line ---- *)
Print["\n=== First ", config["zeroK"], " zeta zeros on the line ==="];
zrc = ZerosOnLineCheck[config["zeroK"]];
Export[FileNameJoin[{exportDir, "zeros_on_line.json"}], zrc, "JSON"];
Print["sign failures: ", zrc["signChangeFailures"], "  up to T = ", N[zrc["tmax"], 10],
  "  elapsed: ", N[zrc["elapsedSeconds"], 4], "s"];

zc = ZeroCountingCheck[config["zeroK"]];
Export[FileNameJoin[{exportDir, "zero_counting.json"}], zc, "JSON"];
Print["N(T) formula vs exact count: rel err ", N[zc["relativeError"], 4]];

figZ = Labeled[
  Plot[RiemannSiegelZ[t], {t, 1, 250}, PlotRange -> All,
   Epilog -> {Red, PointSize[0.012], Point[{#, 0} & /@ N[Im[ZetaZero[Range[18]]]]]},
   AxesLabel -> {"t", "Z(t)"}, PlotLabel -> "Riemann-Siegel Z on the critical line"],
  "First zeros of Z(t): all crossings on the line Re(s)=1/2", Top];
Export[FileNameJoin[{exportDir, "fig_zeros.png"}], figZ, ImageResolution -> 200];

(* ---- 3. Explicit formula: rebuild psi(x) from zeros ---- *)
Print["\n=== Psi reconstruction from ", config["psiK"], " zeros ==="];
psi = PsiReconstruction[config["psiX"], config["psiK"]];
d = psi["data"];
Export[FileNameJoin[{exportDir, "psi_data.csv"}], Prepend[d, {"x", "psi_actual", "psi_reconstructed", "error"}], "CSV"];
figPsi = Labeled[
  ListLinePlot[{d[[All, {1, 2}]], d[[All, {1, 3}]]},
   PlotLegends -> {"psi(x) exact sieve", "reconstructed from " <> ToString[config["psiK"]] <> " zeros"},
   AxesLabel -> {"x", "psi(x)"}, PlotRange -> All],
  "Explicit formula: primes rebuilt from zeta zeros", Top];
Export[FileNameJoin[{exportDir, "fig_psi_reconstruction.png"}], figPsi, ImageResolution -> 200];
figPsiErr = Labeled[
  ListLinePlot[d[[All, {1, 4}]], AxesLabel -> {"x", "error"}, PlotRange -> All],
  "psi error: O(x log^2 x / T), T = " <> ToString[N[psi["maxOrdinate"], 5]], Top];
Export[FileNameJoin[{exportDir, "fig_psi_error.png"}], figPsiErr, ImageResolution -> 200];
Print["psi maxAbsError: ", N[psi["maxAbsError"], 6]];

(* ---- 4. GRH corner: least prime in arithmetic progressions ---- *)
Print["\n=== Least prime in AP, q <= ", config["apQ"], " ==="];
ap = LeastPrimeAPCheck[config["apQ"]];
apExport = KeyDrop[ap, "rows"];
Export[FileNameJoin[{exportDir, "least_prime_ap_summary.json"}], apExport, "JSON"];
rows = ap["rows"];
Export[FileNameJoin[{exportDir, "least_prime_ap.csv"}],
  Prepend[Map[{#q, #a, #leastPrime, #grhBound} &, rows], {"q", "a", "leastPrime", "grhBound"}], "CSV"];
sample = rows[[RandomSample[Range[Length[rows]], UpTo[2000]]]];
figAP = Labeled[
  LogLogPlot[{2 q^2 Log[q]^2}, {q, 2, config["apQ"]},
   Epilog -> {PointSize[0.004], Black, Point[{#q, #leastPrime} & /@ sample]},
   AxesLabel -> {"q", "least prime p(a,q)"}, PlotLegends -> {"GRH bound 2 q^2 log^2 q"},
   PlotLabel -> "Least prime in AP vs GRH bound"],
  "GRH shadow check: every least prime below the bound", Top];
Export[FileNameJoin[{exportDir, "fig_least_prime.png"}], figAP, ImageResolution -> 200];
Print["classes: ", ap["classesChecked"], "  violations: ", ap["grhViolations"],
  "  maxLeastPrime: ", ap["maxLeastPrime"], "  elapsed: ", N[ap["elapsedSeconds"], 4], "s"];

(* ---- 5. Validation report ---- *)
report = StringJoin[
  "# RH-MantiQ validation report\n\n",
  "Scale: ", scale, "\n\n",
  "## Lagarias RH-equivalence\n",
  "- n <= ", ToString[lag["nmax"]], " checked with ", ToString[lag["processorCount"]], " parallel kernels\n",
  "- violations: ", ToString[lag["violations"]], " (0 means the inequality held everywhere; n=1 is the known equality case)\n",
  "- tightest n>1: ", ToString[lag["tightest"][[1]]["n"]], " with ratio ", ToString[N[lag["tightest"][[1]]["ratio"], 8]], "\n",
  "- elapsed: ", ToString[N[lag["elapsedSeconds"], 4]], " s\n\n",
  "## Zeros on the critical line\n",
  "- first ", ToString[zrc["kmax"]], " ordinates checked by sign changes of RiemannSiegelZ between consecutive zeros\n",
  "- sign-change failures: ", ToString[zrc["signChangeFailures"]], "\n",
  "- T max = ", ToString[N[zrc["tmax"], 10]], "\n",
  "- N(T) counting formula relative error: ", ToString[N[zc["relativeError"], 4]], "\n\n",
  "## Explicit formula\n",
  "- psi(x) rebuilt from ", ToString[psi["kzeros"]], " zeros up to x = ", ToString[psi["xmax"]], "\n",
  "- max abs error: ", ToString[N[psi["maxAbsError"], 6]], "\n\n",
  "## GRH: least prime in arithmetic progression\n",
  "- reduced residue classes checked: ", ToString[ap["classesChecked"]], "\n",
  "- GRH bound violations: ", ToString[ap["grhViolations"]], "\n",
  "- largest least prime seen: ", ToString[ap["maxLeastPrime"]], "\n\n",
  "Machine: ", machine["WolframVersion"], ", ", ToString[machine["ProcessorCount"]], " cores\n"];
Export[FileNameJoin[{exportDir, "validation_report.md"}], report, "Text"];
Print["\nBUILD COMPLETE: ", scale, " scale"];
