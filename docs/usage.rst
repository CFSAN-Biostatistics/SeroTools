========
Usage
========

.. highlight:: none

SeroTools provides methods for querying, comparing, and predicting serovar names and antigenic formulae.

.. _query-label:

query
-----

The White-Kauffmann-Le Minor (WKLM) repository can be queried by submitting an appropriate 
serovar name or antigenic formula in an input file composed of a single query per line:: 

    $ serotools query -i <input_file>
    
or as a command line argument::

    $ serotools query -s 'Paratyphi A'
    
Output::

    Input	     Name	      Formula
    Paratyphi A	 Paratyphi A  I [1],2,12:a:[1,5]

.. _compare-label:

compare
-------

SeroTools provides functionality for comparing serovar predictions by evaluating multiple
states of congruence (exact, congruent, minimally congruent, incongruent). Serovar names 
and/or antigenic formulae may be submitted as a tab-delimited input file composed of two 
columns of serovar predictions::  

    $ serotools compare -i <input_file>

or as command line arguments::

    $ serotools compare -1 'Hull' -2 'I 16:b:1,2'

Output::

    Subject	 Name  Formula	   Query	   Name	 Formula	 Result
    Hull	 Hull  I 16:b:1,2  I 16:b:1,2  Hull	 I 16:b:1,2	 exact match

.. _predict-label:

predict
-------

SeroTools also provides the ability to predict which serovars are most likely
represented by (minimally congruent with) a serovar name or antigenic formula, which
may be submitted in an input file composed of a single query on each line::

    $ serotools predict -i <input_file>
    
or as a command line argument::

    $ serotools predict -s 'I 2,12:–:–'

Output::

    Input       MinCon_Name  MinCon_Formula
    I 2,12:–:–  Kiel         I [1],2,12:g,p:–
    I 2,12:–:–  Koessen      I 2,12:l,v:1,5
    I 2,12:–:–  Nitra        I 2,12:g,m:–
    I 2,12:–:–  Paratyphi A  I [1],2,12:a:[1,5]
      
      
Input Formatting

The input data must follow the nomenclature conventions as specified in [1]_.

- Named serovars generally compose a single word (no spaces), with a few exceptions (e.g. Paratyphi B, II Alsterdorf), and should not contain hyphens or additional information. 
- Named variants specified in [1]_ must adhere to the following format: <SerovarName> var. <factor>+,<factor>+. For example, Westhampton var. 15+ or Westhampton var. 15+,34+. 
- In order to query the WKLM database with an antigenic formula (*serotools query*), the formula must match the expected nomenclature including optional '[]', exclusive '{}', or weakly agglutinable '()' factor designations. A space must separate the subspecies designation from the antigens; colons must separate the antigens; commas must separate the antigenic factors. (e.g. I 1,4,[5],12:e,h:1,5:[R1…]). Missing antigens should be specified using '–' (e.g. I 4,12,27:b:– or I 1,9,12:–:–). Input may contain the terms Nonmotile, Rough, or Mucoid, which will be converted into the appropriate format. In addition, phage conversion factors denoted by underlining in [1]_ should be denoted as optional '[]', with the exception of exclusive phage conversion factors (e.g. {15} and {15,34}). 
- Note - *serotools predict* does not require optional factors.


.. [1] Grimont PA, Weill FX. Antigenic Formulae of the Salmonella Serovars. 9th. Paris, France: WHO Collaborating Center for Reference and Research on Salmonella, Institut Pasteur; 2007 <https://www.pasteur.fr/sites/default/files/veng_0.pdf>.

