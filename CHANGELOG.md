# Change Log

## [0.2.1] - 2015-08-06
### Fixed
- Alignment stats shown at completion now properly account for supplementary alignments

## [0.2.0] - 2015-08-02
### Added
- Changelog (history starts here)
- Option to redirect all output (reports and SAM) to file prefix and not stdout
- Print basic alignment stats upon completion
- Show count percentages in UniProt report (both full and basic)

### Fixed
- Option for minimum ORF length temporarily taken out (to be added in version 0.3.0)
- Typo in Uniprot report
- A few minor compilation warnings

### Changed
- Option for minimum ORF length now set to argument "-l".  "-o" used for output file prefix
- Report type "-u" option starts at 0 now, with 1 being default (full report)
- Reorganized usage message to group protein/ORF related options into their own section (to be used more in future development)

### Known Issues
- Alignment stats shown at completion do not properly subtract out supplementary alignments, so total and percentages are slightly off from flagstat output
- File output currently requires prefix to be given, will not autogenerate name if one isn't specified

