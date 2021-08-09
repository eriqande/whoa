## Test environments

* local OS X install, R 4.0.3
* check_win_devel(), check_win_release(), and check_win_oldrelease()
* check_rhub()

## R CMD check results

MAC, local: 0 errors | 0 warnings | 0 notes
WINDOWS, win-builder release: 0 errors | 0 warnings | 1 note
  The 1 note = first submission 
WINDOWS, win-builder devel: 0 errors | 0 warnings | 1 note
  The 1 note = first submission
WINDOWS, win-builder oldrelease: 0 errors | 0 warnings | 1 note
  The 1 note = first submission
WINDOWS, Rhub Windows Server 2008: 0 errors | 0 warnings | 1 notes
  1 note = first submission
LINUX, Rhub Fedora,  0 errors | 0 warnings | 1 note
  1 note = first submission
LINUX, Rhub Ubuntu 20.04,  0 errors | 0 warnings | 1 note
  1 note = first submission


## Downstream dependencies

Currently no known reverse dependencies

## User Notices

* Fixed the datadryad URL, which has been moved.
* This package was previously on CRAN but it was removed because one of its
dependencies (vcfR) had some compile errors that were not fixed in time.
Those compile errors are now fixed with vcfR, and I am resubmitting this to be
back up on CRAN.  Thank you.



