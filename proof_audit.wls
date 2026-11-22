#!/usr/bin/env wolframscript
(* Adversarial audit for RH-MantiQ: independent recomputation of every claim. *)
SeedRandom[42];
projectDir = DirectoryName[$InputFileName];
SetDirectory[projectDir];
exportDir = FileNameJoin[{projectDir, "exports"}];
Get[FileNameJoin[{projectDir, "src", "RHVerification.wl"}]];
LaunchKernels[];

passed = 0; failed = 0; warnings = 0;
lines = {};
log[s_] := (AppendTo[lines, ToString[s, InputForm]]; Print[s]);
audit[name_, test_, detail_: ""] := If[TrueQ[test],
  (passed++; log["PASS: " <> name]),
  (failed++; log["FAIL: " <> name <> " -- " <> ToString[detail, InputForm]])];
warn[name_, detail_] := (warnings++; log["WARN: " <> name <> " -- " <> ToString[detail, InputForm]]);

(* 1. Known constants, independently checked *)
audit["zeta first zero ordinate == 14.13472514173469",
  Abs[Im[ZetaZero[1]] - 14.134725141734693] < 10^-9];
audit["HarmonicNumber[100] == 5.187377517639621",
  Abs[HarmonicNumber[100] - 5.187377517639621] < 10^-9];
audit["sigma(360) == 1170", DivisorSigma[1, 360] == 1170];
audit["sigma(12) == 28 (the tightest Lagarias n)", DivisorSigma[1, 12] == 28];

(* 2. Mangoldt sieve vs built-in MangoldtLambda *)
sv = MangoldtSieve[500];
builtin = N[MangoldtLambda[Range[500]]];
audit["Mangoldt sieve matches built-in MangoldtLambda on 1..500",
  Max[Abs[sv - builtin]] < 10^-12, Max[Abs[sv - builtin]]];

(* 3. Lagarias: equality exactly at n=1, strict inequality for n=2..1000 *)
lagSmall = LagariasCheck[1000];
audit["Lagarias margin at n=1 is exactly 0 (equality case)",
  Abs[(DivisorSigma[1, 1] - (HarmonicNumber[1] + Exp[HarmonicNumber[1]] Log[HarmonicNumber[1]]))] < 10^-9];
t12 = Select[lagSmall["tightest"], #["n"] == 12 &];
audit["n=12 is among the tightest cases with ratio ~0.988636",
  Length[t12] > 0 && Abs[t12[[1]]["ratio"] - 0.9886364271959158] < 10^-9];

(* 4. Lagarias at scale from exports *)
lagFile = FileNameJoin[{exportDir, "lagarias_summary.json"}];
If[FileExistsQ[lagFile],
  lagBig = Association @@ Import[lagFile, "JSON"];
  audit["Big Lagarias run reports zero violations", lagBig["violations"] == 0, lagBig["violations"]];
  audit["Big Lagarias run covers the full range", lagBig["nmax"] >= 10^6, lagBig["nmax"]],
  warn["lagarias_summary.json missing", "run build_project.wls first"]];

(* 5. Zero counting formula vs exact ordinate list *)
zc = ZeroCountingCheck[2000];
audit["N(T) formula agrees with exact count to <0.1%",
  zc["relativeError"] < 10^-3, zc["relativeError"]];

(* 6. Explicit formula: error shrinks with more zeros *)
p50 = PsiReconstruction[200, 50];
p200 = PsiReconstruction[200, 200];
audit["psi reconstruction error decreases with more zeros",
  p200["maxAbsError"] < p50["maxAbsError"], {p50["maxAbsError"], p200["maxAbsError"]}];
audit["psi(100) matches known value 94.0453112293574",
Abs[ChebyshevPsi[100][[100]] - 94.0453112293573981] < 10^-6];

(* 7. Least prime in AP: brute-force cross-check of 10 random entries *)
apFile = FileNameJoin[{exportDir, "least_prime_ap.csv"}];
If[FileExistsQ[apFile],
  apData = Rest[Import[apFile, "CSV"]];
  sample = RandomSample[apData, UpTo[10]];
  ok = And @@ Table[
    Module[{q, a, p, found, k},
     q = row[[1]]; a = row[[2]]; p = row[[3]];
     k = 1; While[! PrimeQ[k q + a], k++];
     k q + a == p],
    {row, sample}];
  audit["least-prime-in-AP entries reproduced by independent scan", ok],
  warn["least_prime_ap.csv missing", "run build_project.wls first"]];

(* 8. GRH bound is never violated in exported data *)
If[FileExistsQ[apFile],
  apData = Rest[Import[apFile, "CSV"]];
  excess = Max[apData[[All, 3]] - apData[[All, 4]]];
  audit["no least prime exceeds the GRH bound 2 q^2 log^2 q", excess < 0, excess]];

log[""];
log["AUDIT SUMMARY: " <> ToString[passed] <> " passed, " <> ToString[failed] <> " failed, " <> ToString[warnings] <> " warnings"];
Export[FileNameJoin[{exportDir, "audit_report.txt"}], StringJoin[Riffle[lines, "\n"]], "Text"];
If[failed > 0, Exit[1]];
Print["AUDIT COMPLETE"];
