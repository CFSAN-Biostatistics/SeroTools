.. :changelog:

History
=======

0.2.1 (2020-09-04)
---------------------

* Updated documentation
* Added JOSS manuscript


0.2.0 (2020-02-17)
---------------------

Significant updates in this version - not backwards compatible.

* The underlying data structures have been converted to pandas Series and DataFrames.
* New 'cluster' subcommand functionality provides the most abundant serovar(s) for clusters of isolates. 
* The 'predict' subcommand functionality has been merged into the 'query' subcommand, such that the default query will return any exact, congruent, and minimally congruent matches unless only exact matches are desired.
* The WKL repository is now available as a pandas DataFrame, in addition to dictionaries and lists.


0.1.1 (2019-11-27)
---------------------

* Corrected a variable name in cli.py
* Updated the algorithm for minimally congruent serovars


0.1.0 (2019-11-19)
---------------------

* Initial version.
