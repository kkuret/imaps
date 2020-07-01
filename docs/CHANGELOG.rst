##########
Change Log
##########

All notable changes to this project are documented in this file.


==================
0.5.0 - 2020-07-01
==================

Changed
-------
- Update kmers.py script:

  - Change the way prtxn are defined
  - Allow using repeat masked genome
  - Improve clustering of top kmers
  - Optimise top motifs clustering


==================
0.4.0 - 2020-05-19
==================

Added
-----
- Add DescriptorSheet object


==================
0.3.0 - 2020-03-29
==================

Added
-----
- Add tagging instructions to contributing docs
- Add Sites analysis

Changed
-------
- Reorganize packet structure
- Change testing tools. No we use black and flake8 instead of pylint.
  Also, all code is now formatted with black


==================
0.2.3 - 2020-03-17
==================

Added
-----
- Add iMaps assests values in constants folder
- Add check_assets script

Changed
-------
- Perform subsampling in kmers script if necessary (lowering memory footprint)
- Update naming of output files in kmer script


==================
0.2.2 - 2019-12-16
==================

Added
-----
- Add contributing guidelines for proposing a change

Changed
-------
- Shift distributions by -1 to make kmer positions more intuitive


==================
0.2.1 - 2019-12-11
==================

Fixed
-----
- Another set of smaller fixes in kmers analysis


==================
0.2.0 - 2019-12-09
==================

Added
-----
- Add iCount peak coverage analysis to sandbox

Fixed
-----
- Multiple smaller fixes in kmers analysis


==================
0.1.0 - 2019-12-04
==================

Added
-----
- Add kmers analysis to imaps/sandbox
- Add contributing docs
