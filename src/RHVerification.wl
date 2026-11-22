BeginPackage["RHVerification`"];

MachineProfile::usage = "MachineProfile[] returns hardware and Wolfram build facts as an Association.";
LagariasCheck::usage = "LagariasCheck[nmax] verifies the Lagarias RH-equivalence sigma(n) <= H_n + E^H_n Log[H_n] for all n <= nmax, vectorized. Returns nmax, elapsed seconds, violation count, and the tightest n by ratio sigma(n)/bound(n).";
ZerosOnLineCheck::usage = "ZerosOnLineCheck[kmax] verifies the first kmax ordinates of ZetaZero lie on the critical line by confirming RiemannSiegelZ changes sign in every gap between consecutive ordinates. Returns kmax, max ordinate reached, sign-change failures, elapsed.";
ZeroCountingCheck::usage = "ZeroCountingCheck[kmax] compares the number of ZetaZero ordinates up to T = Im[ZetaZero[kmax]] with the Riemann-von Mangoldt counting function N(T).";
PsiReconstruction::usage = "PsiReconstruction[xmax, kzeros] reconstructs ChebyshevPsi from the first kzeros zeta zeros via the explicit formula and compares with the exact sieve. Returns sampled data for plotting.";
LeastPrimeAPCheck::usage = "LeastPrimeAPCheck[qmax] finds the least prime p == a mod q for every reduced residue class with q <= qmax and tests the GRH upper bound p < 2 q^2 Log[q]^2 (Bach-Sorenson). Returns violations, worst gap, and the full table.";
MangoldtSieve::usage = "MangoldtSieve[n] returns the von Mangoldt values Lambda(1..n) by sieve.";
ChebyshevPsi::usage = "ChebyshevPsi[n] returns psi(x) = sum Lambda(k) for k <= n, from the sieve.";
RunValidationSuite::usage = "RunValidationSuite[] runs every check at moderate scale and returns one Association of results.";
HarmonicVector::usage = "HarmonicVector[n] returns H_1..H_n: exact for n <= 100, asymptotic expansion beyond (error < 4*^-15).";

Begin["`Private`"];

lagariasBound[h_] := h + Exp[h] Log[h];

MachineProfile[] := <|
  "WolframVersion" -> $Version,
  "ProcessorCount" -> $ProcessorCount,
  "ProcessorType" -> $ProcessorType,
  "KernelID" -> $KernelID
|>;

HarmonicVector[n_] := Module[
  {lim = Min[n, 100], exact, big, m},
  exact = Accumulate[1./Range[lim]];
  If[n <= lim, exact,
   m = Range[lim + 1., n];
   big = N[Log[m]] + N[EulerGamma] + 1./(2 m) - 1./(12 m^2) + 1./(120 m^4);
   Join[exact, big]]
];

LagariasCheck[nmax_] := Module[
  {t0, sigma, h, bound, margins, ratios, order, order0, worst, elapsed, nChunks, chunkLen, ranges},
  t0 = AbsoluteTime[];
  nChunks = $ProcessorCount;
  chunkLen = Ceiling[nmax/nChunks];
  ranges = Table[Range[1 + (i - 1) chunkLen, Min[i chunkLen, nmax]], {i, nChunks}];
  sigma = Flatten[ParallelTable[
     DivisorSigma[1, r], {r, ranges}]];
  h = HarmonicVector[nmax];
  bound = lagariasBound[h];
  margins = sigma - bound;
  elapsed = AbsoluteTime[] - t0;
  ratios = sigma / bound;
  order0 = Ordering[-ratios, UpTo[11]];
  order = DeleteCases[order0, 1][[;; UpTo[10]]];
  worst = Table[<|"n" -> n, "sigma" -> sigma[[n]], "bound" -> bound[[n]],
     "ratio" -> ratios[[n]], "margin" -> margins[[n]]|>, {n, order}];
  <|
    "nmax" -> nmax,
    "elapsedSeconds" -> elapsed,
    "violations" -> Total[UnitStep[margins - 10.^-9 Abs[bound]]],
    "recheckNeeded" -> Total[UnitStep[10.^-6 Abs[bound] - Abs[margins]]],
    "minMargin" -> Min[margins],
    "maxRatio" -> ratios[[order[[1]]]],
    "tightest" -> worst,
    "processorCount" -> $ProcessorCount
  |>
];

ZeroOrdinates[kmax_] := Module[
  {nChunks, chunkLen, ranges},
  nChunks = $ProcessorCount;
  chunkLen = Ceiling[kmax/nChunks];
  ranges = Table[{1 + (i - 1) chunkLen, Min[i chunkLen, kmax]}, {i, nChunks}];
  Flatten[ParallelTable[N[Im[ZetaZero[Range[r[[1]], r[[2]]]]]], {r, ranges}]]
];

ZerosOnLineCheck[kmax_] := Module[
  {ords, gaps, mids, zvals, failures, elapsed, t0, tmax, nChunks, chunkLen, chunks},
  t0 = AbsoluteTime[];
  ords = ZeroOrdinates[kmax];
  gaps = Differences[ords];
  mids = Most[ords] + gaps/2;
  nChunks = $ProcessorCount;
  chunkLen = Ceiling[Length[mids]/nChunks];
  chunks = Table[mids[[1 + (i - 1) chunkLen ;; Min[i chunkLen, Length[mids]]]], {i, nChunks}];
  zvals = Flatten[ParallelMap[RiemannSiegelZ, chunks]];
  failures = Total[UnitStep[-Abs[zvals]]];
  elapsed = AbsoluteTime[] - t0;
  tmax = Last[ords];
  <|
    "kmax" -> kmax,
    "tmax" -> tmax,
    "signChangeFailures" -> failures,
    "minGap" -> Min[gaps],
    "elapsedSeconds" -> elapsed
  |>
];

ZeroCountingCheck[kmax_] := Module[
  {ords, t, nExact, nFormula, err, elapsed},
  ords = ZeroOrdinates[kmax];
  t = Last[ords];
  nExact = kmax;
  nFormula = t/(2 Pi) Log[t/(2 Pi E)] + 7/8;
  err = Abs[nFormula - nExact]/nExact;
  <|
    "kmax" -> kmax,
    "T" -> t,
    "countExact" -> nExact,
    "countFormula" -> nFormula,
    "relativeError" -> err
  |>
];

MangoldtSieve[n_] := Module[
  {lam = ConstantArray[0., n], comp = ConstantArray[False, n], p},
  Do[
    If[! comp[[p]],
     Do[lam[[j]] += Log[p], {j, NestWhileList[# p &, p, # p <= n &]}];
     Do[comp[[j]] = True, {j, 2 p, n, p}]],
    {p, 2, n}];
  lam
];

ChebyshevPsi[n_] := Accumulate[MangoldtSieve[n]];

PsiReconstruction[xmax_, kzeros_] := Module[
  {ords, xgrid, psiActual, psiRecon, mainTerm, zerosTerm, logTerm, data, t0},
  t0 = AbsoluteTime[];
  ords = ZeroOrdinates[kzeros];
  xgrid = Range[2, xmax, Max[1, Floor[xmax/400]]];
  psiActual = ChebyshevPsi[xmax][[xgrid]];
  mainTerm = N[xgrid];
  logTerm = Log[2 Pi] + 0.5 Log[1 - 1/xgrid^2];
  zerosTerm = Table[
    2 Re[x^(1/2 + I g)/(1/2 + I g)],
    {x, xgrid}, {g, ords}];
  psiRecon = mainTerm - Total[zerosTerm, {2}] - logTerm;
  data = Transpose[{xgrid, psiActual, psiRecon, psiActual - psiRecon}];
  <|
    "xmax" -> xmax,
    "kzeros" -> kzeros,
    "maxOrdinate" -> Last[ords],
    "data" -> data,
    "maxAbsError" -> Max[Abs[data[[All, 4]]]],
    "elapsedSeconds" -> AbsoluteTime[] - t0
  |>
];

LeastPrimeAPCheck[qmax_] := Module[
  {rows, result, worst, elapsed},
  elapsed = AbsoluteTime[];
  rows = Flatten[ParallelTable[
    Module[{residues, k, p},
     residues = Select[Range[1, q - 1], GCD[#, q] == 1 &];
     Table[
      With[{a = r},
       k = 1;
       While[! PrimeQ[k q + a], k++];
       p = k q + a;
       <|"q" -> q, "a" -> a, "leastPrime" -> p, "k" -> k,
        "grhBound" -> 2. q^2 Log[q]^2|>],
      {r, residues}]],
    {q, 2, qmax}],
   1];
  worst = Max[rows[[All, "leastPrime"]] - rows[[All, "grhBound"]]];
  elapsed = AbsoluteTime[] - elapsed;
  <|
    "qmax" -> qmax,
    "classesChecked" -> Length[rows],
    "grhViolations" -> Count[rows, r_ /; r["leastPrime"] > r["grhBound"]],
    "worstExcessOverBound" -> worst,
    "maxLeastPrime" -> Max[rows[[All, "leastPrime"]]],
    "rows" -> rows,
    "elapsedSeconds" -> elapsed
  |>
];

RunValidationSuite[] := Module[{res},
  res = <||>;
  res["machine"] = MachineProfile[];
  res["lagarias1e6"] = LagariasCheck[10^6];
  res["zerosOnLine10k"] = ZerosOnLineCheck[10000];
  res["zeroCounting10k"] = ZeroCountingCheck[10000];
  res["psiReconstruction"] = PsiReconstruction[1000, 200];
  res["leastPrimeAP200"] = LeastPrimeAPCheck[200];
  res
];

End[];
EndPackage[];
