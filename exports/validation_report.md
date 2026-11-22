# RH-MantiQ validation report

Scale: full

## Lagarias RH-equivalence
- n <= 100000000 checked with 24 parallel kernels
- violations: 0 (0 means the inequality held everywhere; n=1 is the known equality case)
- tightest n>1: 12 with ratio 0.988636
- elapsed: 51.37 s

## Zeros on the critical line
- first 20000 ordinates checked by sign changes of RiemannSiegelZ between consecutive zeros
- sign-change failures: 0
- T max = 18046.5
- N(T) counting formula relative error: 0.0000300927

## Explicit formula
- psi(x) rebuilt from 1000 zeros up to x = 2000
- max abs error: 4.72932

## GRH: least prime in arithmetic progression
- reduced residue classes checked: 304191
- GRH bound violations: 0
- largest least prime seen: 70783

Machine: 15.0.1 for Microsoft Windows (64-bit) (July 2, 2026), 24 cores
