## Test environments

* local OS X install, R 3.5.1, 3.6.0
* ubuntu 14.04 (on travis-ci), R 3.5.0
* check_win_devel() and check_win_release()
* check_rhub()

## R CMD check results

MAC, local: 0 errors | 0 warnings | 0 notes
LINUX, travis-ci: 0 errors | 0 warnings | 0 notes
WINDOWS, win-builder release: 0 errors | 0 warnings | 1 note
  The 1 note = first submission
WINDOWS, win-builder devel: 0 errors | 0 warnings | 1 note
  The 1 note = first submission
WINDOWS, Rhub Windows Server 2008: 0 errors | 0 warnings | 2 notes
  1 note = first submission
  1 note = the example output directories were left in the package.  This seems to be
  an inconsistency on the Rhub side of things?
LINUX, Rhub Fedora,  0 errors | 0 warnings | 1 note
  1 note = first submission
LINUX, Rhub Ubuntu 16.04,  0 errors | 0 warnings | 1 note
  1 note = first submission
LINUX, Rhub Debian,  0 errors | 0 warnings | 0 note

## Downstream dependencies

Currently no known reverse dependencies

## User Notices

This is the first submission of this package.  



