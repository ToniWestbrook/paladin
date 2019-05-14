# Change Log

## [1.4.6] - 2019-05-07
### Fixed
- Web request retries now work properly

### Added
- Taxonomic ID and Lineage to UniProt report

## [1.4.5] - 2019-05-07
### Changed
- Changed UniProt web request batch size to 5000 entries to reduce chances of server timeout

## [1.4.4] - 2018-09-12
### Fixed
- Corrected issues causing compiler warnings

## [1.4.3] - 2018-07-25
### Changed
- Changed all web requests from HTTP to HTTPS to accommodate UniProt's new requirements

## [1.4.2] - 2018-03-20
### Fixed
- Corrected issue that sometimes caused PALADIN to hang when downloading data from UniProt
- Corrected ability to recover from UniProt server errors

## [1.4.1] - 2017-04-17
### Fixed
- Corrected truncated translation code descriptions

## [1.4.0] - 2017-03-23
### Added
- Option to detect ORFs/translate/align across multiple non-standard genetic codes (-z option). Details of available translation codes may be found here: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

## [1.3.2] - 2017-02-07
### Added
- Alignment command can now directly take a protein multi-FASTA and skip ORF detection (-p option)
- Prepare and alignment commands can now make use of a proxy server (HTTP/SOCKS) for contacting UniProt (-P option)

## [1.3.1] - 2017-01-01
### Added
- UniProt report now includes max mapping quality for each protein

### Changed
- UniProt report floating precision now set to 5 digits

## [1.3.0] - 2016-10-16
### Added
- UniProt report now includes average mapping quality for each protein

### Fixed
- Corrected mapping quality calculations to better reflect probability in protein space
- Bug that would arise during report generation when UniProt servers were busy, will retry now

### Changed
- Default alignment threshold score to 15, mismatch penalty to 3, clipping penalty to 0, and open gap penalty to 0 (in accordance with empirical testing)

## [1.2.1] - 2016-04-05
### Changed
- All file IO now checks for open errors (a few weren't checking previously)

## [1.2.0] - 2016-03-27
### Added
- New columns in UniProt report related to database cross-references (KEGG, NCBI, PATRIC, Ensembl)
- Added experimental support/fix for shared memory indexed references

## [1.1.0] - 2015-11-08
### Added
- Option for preparing a pre-downloaded and/or pre-indexed reference.  This will skip the download and clean any applicable portion of the reference/index
- Index version compatibility system.  Will ensure at alignment time the index is compatible with current software version
- Additional runtime status reporting, especially in areas where processing can take a significant amount of time (ORF detection, index preparation, index loading, etc)

### Changed
- UniRef50 option to UniRef90 (in accordance with empirical testing)
- A few minor description changes

## [1.0.3] - 2015-09-30
### Added
- UniRef50 as one of the preparation/download options

## [1.0.2] - 2015-09-29
### Added
- UniProt reporting now properly parses header styles of both SwissProt/TrEMBL and UniRef

## [1.0.1] - 2015-09-28
### Changed
- Default alignment threshold score from 30 to 20 (in accordance with empirical testing)

## [1.0.0] - 2015-08-11
### Added
- Minimum ORF length filtering
- Related options (Constant minimum length, percentage minimum length, adjustment for smaller read lengths)

### Fixed
- Translating edges of detected ORFs (previously were truncated)

### Changed
- A few minor description changes

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

