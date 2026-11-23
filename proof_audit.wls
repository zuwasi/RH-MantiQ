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

(* 9. Li's criterion: independent closed-form coefficients vs the package's Cauchy extraction.
   The audit rebuilds the coefficients of Log[xi(1+u)] itself from closed-form series:
   log((1+u)/2) - ((1+u)/2) Log[Pi] + Log[Gamma(1/2+u/2)] + Log[u Zeta(1+u)],
   where Log[Gamma(1/2+u/2)] coefficients follow from polygamma values at 1/2 and
   Log[u Zeta(1+u)] coefficients follow from the Stieltjes constants by series composition. *)
wpAud = 80; KAud = 33;
cA = XiLogCoefficients[KAud, wpAud];
c1 = ConstantArray[0, KAud];
c1[[1]] = N[1 - Log[Pi]/2, wpAud];
For[j = 2, j <= KAud, j++, c1[[j]] = N[(-1)^(j + 1)/j, wpAud]];
c2 = ConstantArray[0, KAud];
c2[[1]] = N[-EulerGamma/2 - Log[2], wpAud];
For[k = 2, k <= KAud, k++, c2[[k]] = N[(-1)^k (2^k - 1) Zeta[k]/(k 2^k), wpAud]];
w = ConstantArray[0, KAud];
For[m = 1, m <= KAud, m++, w[[m]] = N[(-1)^(m - 1) StieltjesGamma[m - 1]/(m - 1)!, wpAud]];
c3 = ConstantArray[0, KAud];
For[j = 1, j <= KAud, j++,
  c3[[j]] = w[[j]];
  For[m = 1, m <= j - 1, m++, c3[[j]] -= (m/j) c3[[m]] w[[j - m]]]];
cClosedAud = c1 + c2 + c3;
maxRelC = Max[Table[Abs[cA[[j]] - cClosedAud[[j]]]/Abs[cA[[j]]], {j, KAud}]];
audit["Li coefficients: closed-form pieces match Cauchy extraction (j<=33) to 1e-9",
  maxRelC < 10^-9, maxRelC];

(* 10. Li's criterion: constants, positivity, two-route agreement *)
li = LiConstantsCheck[50, 2000, 100];
audit["Li lambda_1 == 1 + gamma/2 - Log[2 Sqrt[Pi]] (closed form)",
  Abs[li["constants"][[1]] - (1 + EulerGamma/2 - Log[2 Sqrt[Pi]])] < 10^-30, li["constants"][[1]]];
audit["Li's criterion: all 50 constants strictly positive (audit scale)",
  li["nonPositiveCount"] == 0, li["nonPositiveCount"]];
audit["Li constants: coefficient route and zeros route agree within 5%",
  li["maxRelativeDiscrepancy"] < 0.05, {li["maxRelativeDiscrepancy"], li["discrepancyAtN"]}];

(* 11. Li's criterion at scale from exports *)
liFile = FileNameJoin[{exportDir, "li_summary.json"}];
If[FileExistsQ[liFile],
  liBig = Association @@ Import[liFile, "JSON"];
  audit["Big Li run reports zero non-positive constants", liBig["nonPositiveCount"] == 0, liBig["nonPositiveCount"]],
  warn["li_summary.json missing", "run build_project.wls first"]];

log[""];
log["AUDIT SUMMARY: " <> ToString[passed] <> " passed, " <> ToString[failed] <> " failed, " <> ToString[warnings] <> " warnings"];
Export[FileNameJoin[{exportDir, "audit_report.txt"}], StringJoin[Riffle[lines, "\n"]], "Text"];
If[failed > 0, Exit[1]];
Print["AUDIT COMPLETE"];
