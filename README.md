# Structural Multiplicity Without Dimensional Proliferation

**Density, Admissibility, and Four-Dimensional Sufficiency in TCGS–SEQUENTION**

Henry Arellano-Peña · Manuscript dated 12 October 2026

This repository contains the numerical verification code, executed results, and supporting audit documents accompanying the final 48-page mathematical and methodological manuscript. No new experiment was performed, and the synthetic calculations are not claimed as empirical validation of Counterspace.

## Deposit status

The verification materials are deposited. The final PDF is prepared but has not yet been transferred to this repository. The remaining file is `main.pdf` (613,423 bytes), with SHA-256:

```text
b5829446a71f1a2c30f82da2235e38c343732f1d90f894c8a774635e5776fe49
```

Add that exact PDF at the repository root to complete the deposit. The numerical verifier does not depend on the PDF.

LaTeX manuscript sources are intentionally excluded, along with source-containing patches, Overleaf archives, and historical source directories. The `.gitignore` reinforces this exclusion for future ordinary Git additions.

## Reproducibility materials

- [Complete verification script](verify_examples.py)
- [Executed verification results](verified_examples.json)
- [Coverage and execution report](verification_report.md)
- [Targeted reference verification](reference_verification.json)
- [Reference audit](reference_audit.md)
- [Original build and preservation audit](validation_report.json)
- [Checksums for the deposited verification materials](SHA256SUMS.txt)

## Run all verification checks

Use Python 3.10 or newer and the standard library:

```sh
python verify_examples.py
```

Do not enable `-O` or `PYTHONOPTIMIZE`, because the verifier uses assertions. The script runs the original examples and all three amendment suites. Its entry point invokes `verify_amendments`, `verify_second_amendment`, and `verify_third_amendment`. Results are written to `verified_examples.json` beside the script. No network access, external dataset, or third-party numerical package is required.

The complete release script has SHA-256:

```text
c6bafd9a4ac4fa0948c61a402a03c05a3dc41344d2658033d27b440d53570a5a
```

The checks cover the printed synthetic tables, exact graph responses, the common-kernel construction, loss-aware accessibility calibration, generic block identities, unequal shared splitting, scale blindness, accessibility-only closure, uncertainty bounds, and separating instruments. Finite numerical checks supplement the manuscript's proofs; they do not replace them. Minor last-bit floating-point differences can occur between platforms.

## Interpretation and scope

Zero-phase relative calibration tests the broader identical-shared/orthogonal-private block structure. A phase scan tests equal shared splitting within that class. Independently corrected absolute accessibility tests the normalization consequence. The accessibility-only closure is an additional diagnostic, not statistically independent confirmation obtained by reusing the same counts. No unique universal ECL or physical four-dimensional source is inferred from a passing finite test.

The verification and validation reports record the earlier local release audit. Statements in those historical reports that no remote deposit had yet occurred refer to that audit session. The present deposit is established by this repository's actual commit history; the PDF remains pending as stated above.

## License

Copyright © 2026 Henry Arellano-Peña. The authored manuscript and accompanying materials are distributed under the Creative Commons Attribution 4.0 International License (CC BY 4.0): https://creativecommons.org/licenses/by/4.0/ . External scientific sources retain their respective attribution requirements.
