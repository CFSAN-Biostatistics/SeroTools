=================
Input Formatting
=================

.. highlight:: none

 
The input data must follow the nomenclature conventions as specified in [1]_.

- Named serovars generally compose a single word (no spaces), with a few exceptions (e.g. Paratyphi B, II Alsterdorf), and should not contain hyphens or additional information. 
- Named variants specified in [1]_ must adhere to the following format: <SerovarName> var. <factor>+,<factor>+. For example, Westhampton var. 15+ or Westhampton var. 15+,34+. 
- In order to receive exact matches when querying the WKLM database with an antigenic formula (*serotools query*), the formula must match the expected nomenclature including optional '[]', exclusive '{}', or weakly agglutinable '()' factor designations. A space must separate the subspecies designation from the antigens; colons must separate the antigens; commas must separate the antigenic factors. (e.g. I 1,4,[5],12:e,h:1,5:[R1…]). Missing antigens should be specified using '–' (e.g. I 4,12,27:b:– or I 1,9,12:–:–). Input may contain the terms Nonmotile, Rough, or Mucoid, which will be converted into the appropriate format. In addition, phage conversion factors denoted by underlining in [1]_ should be denoted as optional '[]', with the exception of exclusive phage conversion factors (e.g. {15} and {15,34}). 

.. [1] Grimont PA, Weill FX. Antigenic Formulae of the Salmonella Serovars. 9th. Paris, France: WHO Collaborating Center for Reference and Research on Salmonella, Institut Pasteur; 2007 <https://www.pasteur.fr/sites/default/files/veng_0.pdf>.
