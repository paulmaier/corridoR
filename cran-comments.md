## Resubmission

This is a resubmission. In this version I have:

* Put software and file-format names in the Description in single quotes,
  and replaced "FST" with "fixation indices".
* Wrapped the `corridor_extract()` example, which took more than 5 s on
  Debian, in `\donttest{}`.

The remaining spelling NOTE flags author names in references ("Maier et
al.", "Sundqvist et al."); these are correct.

## Test environments

* local macOS (aarch64), R 4.5.3
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)
* win-builder and Debian incoming pre-tests (previous submission)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.
