===============================
SeroTools
===============================


.. Image showing the PyPI version badge - links to PyPI
.. image:: https://img.shields.io/pypi/v/serotools.svg
        :target: https://pypi.python.org/pypi/serotools

.. Image showing the Travis Continuous Integration test status, commented out for now
.. .. image:: https://img.shields.io/travis/CFSAN-Biostatistics/serotools.svg
..        :target: https://travis-ci.org/CFSAN-Biostatistics/serotools

This package serves as a toolkit and repository for the White-Kauffmann-Le Minor scheme for *Salmonella* serotyping, which is made available in multiple formats, along with methods for querying and comparing serovar names and antigenic formulae, as well as determining the most abundant serovar for a cluster of isolates.

SeroTools was developed by the United States Food and Drug Administration, Center for Food 
Safety and Applied Nutrition.

* Free software
* Documentation: https://serotools.readthedocs.io
* Source Code: https://github.com/CFSAN-Biostatistics/serotools
* PyPI Distribution: https://pypi.python.org/pypi/serotools

Introduction
------------

*Salmonella* bacteria are major foodborne pathogens estimated by the U.S. Centers for Disease Control and Prevention to cause 1.35 million infections annually in the United States [1]_. Serological subtyping (serotyping) of *Salmonella* has historically been a critical component of characterization and successful outbreak identification and traceback efforts employed by public health researchers and regulatory agencies. The White-Kauffmann-Le Minor (WKL) *Salmonella* serotyping scheme specifies the commonly accepted naming and formatting conventions for *Salmonella* serotyping data and the antigenic factors (and other characteristics) which define each serovar. *Salmonella* serotyping data is routinely employed by a broad range of scientific researchers, physicians, public health professionals, food safety experts, etc.

SeroTools addresses multiple critical needs for the efficient analysis of *Salmonella* serotyping data. As technological advances continue to produce a range of high resolution subtyping options, including *in silico* serovar prediction based on whole genome sequencing, new tools are necessary for efficient method-comparison studies and quality control applied to increasingly large numbers of isolates. SeroTools serves as the only multiformat WKL repository accessible for software development and provides the only existing tools for querying the WKL scheme, comparing serovars for congruence, and predicting the most abundant serovar for clusters of isolates.

.. [1] The U.S. Centers for Disease Control and Prevention. <https://www.cdc.gov/salmonella/index.html>.


Features
--------

* Query the White-Kauffmann-Le Minor *Salmonella* serotyping repository

* Compare serovar predictions for state of congruence

* Determine the most abundant serovar for a cluster of isolates


Citing SeroTools
--------------------------------------

To cite SeroTools, please reference the SeroTools GitHub repository:

    https://github.com/CFSAN-Biostatistics/serotools


License
-------

See the LICENSE file included in the SeroTools distribution.




