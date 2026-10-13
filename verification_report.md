# Verification script: final release audit

## Executed artifact

The root `verify_examples.py` is **byte-identical** to the complete third-amendment script. It is not the shorter baseline script. It was inspected in full and copied to an otherwise empty directory for a fresh standalone run under Python 3.13.5, with optimization disabled (`sys.flags.optimize == 0`). Do not use `python -O`: this verifier uses assertions.

- Script SHA-256: `c6bafd9a4ac4fa0948c61a402a03c05a3dc41344d2658033d27b440d53570a5a`
- Result SHA-256: `328ebdc87da169cc0b04da43d4f48e9107c59d19cc1dffa2f3938e6b8dc0489c`
- The complete JSON output is equal to the third-amendment output in every key and value.
- The standalone stdout JSON and the written `verified_examples.json` agree.
- The top-level `__main__` entry point calls `main()`, which invokes `verify_amendments(graphs)`, `verify_second_amendment()`, and `verify_third_amendment()` after the original examples.

## Coverage of the printed tables

| Manuscript target | Implementation / JSON location | Final audit |
|---|---|---|
| Table 2: generalized dimensions | `generalized_dimension`; `multifractal_example` | All printed D1 and D2 values checked at six decimals; D0 is four. |
| Table 3: graph responses | `laplacian`, exact `determinant` and `solve`; `graphs` | Both characteristic polynomials and each rational response checked. |
| Table 4: torus spectral dimension | `theta_f`, `spectral_dimension`; `spectral_example` | Every printed six-decimal value checked; cutoff 45 versus 90 checked. |
| Table 5: finite joint kernel | `verify_amendments`; `amendment.joint_kernel_exact_rows` | All seven rows present; exact complex quarter-turn components and squared visibility checked. |

The original four-ball comparison, graph triangles, source matrices and response gaps remain in the output. Finite enumeration is not substituted for the article's integer-form proof.

## All three amendment suites are included and executed

| Suite | Selected completed checks | Recorded maximum floating-point residual |
|---|---|---|
| First amendment | 72 exact quarter-turn cases; 357 complex kernel cases; 7,140 direct port-probability cases; exact label-free resolvent traces; normalized failure completion and held-out test | 6.661338147750939e-16 |
| Second amendment | 136 common-efficiency cases; 544 calibrated-mismatch cases; 680 interval-containment cases; 68 atomic-ball cases; 1,768 loss-completed instruments; 1,001 stability comparisons; invalid-input and domain checks | 3.3306690738754696e-16 |
| Third amendment | 63 exact generic zero-phase cases; 1,764 exact shared-split phase cases; 3,087 complex-vector cases; 11,466 complex and 234 exact rescaling cases; 63 closure-interval cases; 336 corrected-count closures; three separating instruments | 4.460365043975677e-16 |

These are synthetic tests, not statistical experimental evidence. Rational identities use `fractions.Fraction`; general phase/isometry checks use standard-library complex arithmetic. Finite grids supplement the manuscript's proofs and do not establish a theorem over an uncountable parameter domain by sampling.

The third suite contains explicit metadata correcting the attribution of the older zero-phase rejection: the old numerical result is retained, but it belongs already to the weaker identical-shared/orthogonal-private block class. A passing zero-phase ratio test is not reported as unique kernel identification. The absolute closure uses independently corrected accessibility rather than inserting the kernel inverse on both sides, and reused measurements are not declared statistically independent.

## Reproduction

From the extracted release directory:

```sh
python verify_examples.py
```

This overwrites `verified_examples.json` beside the script. Python 3.10 or newer and its standard library are sufficient; no network, external dataset, numerical library, or historical file is loaded. Minor last-bit differences can occur between platforms; exact rational checks and printed table precision remain the intended checks. The supplied checksums identify the exact artifacts used in this release, not a guarantee of cross-platform byte-identical floating-point output.

## What was not done

No experimental anchor was fitted or reanalyzed. No new proof was added in this editorial pass. No remote repository upload, commit, release, journal submission or preprint posting was performed. The repository-deposit package is prepared locally and must not be described as already published.
