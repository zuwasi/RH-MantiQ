# RH-MantiQ: Wolfram Language verification lab for the Riemann Hypothesis news

On November 20, 2026, Anthropic mathematician Levent Alpöge announced a proof of the
Riemann Hypothesis (and the Generalized Riemann Hypothesis), produced with Claude and
released only as formalized Lean code: 14 archive parts, on the order of 17 million lines
(https://alpo.ge/remark.pdf, https://github.com/CocoSgt/riemann-hypothesis-lean).

This project is the Wolfram Language answer to that news: you cannot re-prove the
Riemann Hypothesis in Mathematica, but you can put the hypothesis and its famous
equivalences and consequences on a numerical test bench, and push the checks far
beyond what a casual session does, using all cores of the machine.

## What it checks

1. **Lagarias criterion** (an RH equivalence, Lagarias 2002):
   `sigma(n) <= H_n + E^H_n Log[H_n]` for every n >= 1, where `sigma` is the divisor sum
   and `H_n` the harmonic number. Equality holds only at n=1. We check it vectorized and
   parallelized over all cores up to n = 10^8, and report the tightest cases
   (they are the highly composite numbers: 12, 120, 60, 2520, 5040, ...).
2. **Zeros on the critical line**: the first k zeta zero ordinates from `ZetaZero`, each
   confirmed by a sign change of `RiemannSiegelZ` in the gap to its neighbor, plus the
   Riemann-von Mangoldt counting function `N(T)` cross-check.
3. **Explicit formula**: rebuild `psi(x)` (Chebyshev psi) from the main term plus the
   first k zeta zeros via `psi(x) = x - sum_rho x^rho/rho - Log[2 Pi] - Log[1-1/x^2]/2`,
   and compare against an independent von Mangoldt sieve. This is the concrete
   prime-zeros connection the proof exploits.
4. **GRH corner**: least prime `p == a mod q` for every reduced residue class with
   q <= qmax, checked against the Bach-Sorenson GRH upper bound `p < 2 q^2 Log[q]^2`.
5. **Li's criterion**: the constants
   `lambda_n = Sum over zeta zeros rho of (1 - (1 - 1/rho)^n)`, n = 1..kmax, are
   positive for every n if and only if the Riemann Hypothesis holds. Each constant is
   computed two independent ways: from the coefficients c_j of `Log[xi(1+u)]`
   (extracted by a spectrally exact Cauchy integral on |u|=1 at working precision up
   to 140 digits, then transformed by `lambda_n = n Sum_j C(n-1,j-1) c_j`), and
   directly from the first k zero ordinates as a conditionally convergent paired sum. The two routes
   must agree, and every constant must be strictly positive.

## Files

- `src/RHVerification.wl` - the package (public API: `MachineProfile`, `LagariasCheck`,
  `ZerosOnLineCheck`, `ZeroCountingCheck`, `PsiReconstruction`, `LeastPrimeAPCheck`,
  `LiConstantsCheck`, `XiLogCoefficients`, `MangoldtSieve`, `ChebyshevPsi`,
  `RunValidationSuite`)
- `build_project.wls` - runs all checks and writes every export. Scale via environment
  variable `RH_SCALE` (`demo` = minutes-scale, `full` = 10^8 / 20000 zeros / q<=1000 /
  200 Li constants)
- `proof_audit.wls` - adversarial audit: known constants, sieve vs built-in
  `MangoldtLambda`, equality case, N(T) accuracy, error monotonicity, brute-force
  cross-check of least-prime entries, Li coefficients rebuilt independently from
  closed-form series and cross-checked against the package's Cauchy extraction
- `exports/` - JSON summaries, CSV tables, PNG figures, `validation_report.md`,
  `audit_report.txt`, `build_log.txt`

## Running

```powershell
$env:RH_SCALE = "full"
& "C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe" -file .\build_project.wls
& "C:\Program Files\Wolfram Research\WolframScript\wolframscript.exe" -file .\proof_audit.wls
```

## Status

On November 20, 2026, Levent Alpöge announced the proof; on November 21, 2026, the
Lean and Mathlib communities certified it correct in a joint statement
(https://leanprover-community.github.io/blog/posts/riemann-hypothesis/). This project
is the independent, human-scale numerical view of the same arithmetic: the five
checks below were stress-tested on one workstation while the formal verdict
landed, and both worlds agree.

Demo video (3 minutes, all chapters timestamped):
https://www.youtube.com/watch?v=WwcqxcO8xRU

## Honest scope

- Numerical checks are evidence, not proof. They run on machine precision (margins are
  orders of magnitude larger than rounding error; the audit counts any margin within
  10^-6 relative of the bound as "recheck needed").
- The Lean formalization is the actual proof instrument; this project is the
  independent, human-scale numerical cross-examination.

## Machine

Wolfram 15.0.1, Intel Core Ultra 9 275HX, 24 cores, 63.4 GB RAM.

## License

MIT. See `LICENSE`. The repository intentionally contains only Wolfram Language
source; it includes none of the Lean 4 or mathlib files from the claimed proof.
