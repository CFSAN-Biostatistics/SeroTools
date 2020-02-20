=================
Input Formatting
=================

.. highlight:: none

 
The input data must follow the nomenclature conventions as specified in [1]_.

- Named serovars (subsp. *enterica*) generally compose a single word or concatenation (e.g. Saintpaul) with no whitespace, with a few exceptions (e.g. Paratyphi B, II Alsterdorf), and should not contain hyphens or additional information (e.g. Choleraesuis is correct, while Cholerae-suis and Cholerae suis are incorrect).

- Named variants specified in [1]_ must adhere to the format below. For example, Westhampton var. 15+ or Westhampton var. 15+,34+:: 

    <SerovarName> var. <factor>+,<factor>+

- The proper convention for an antigenic formula is as follows: a space must separate the subspecies symbol from the antigens; colons must separate the antigens; commas must separate the antigenic factors (e.g. I 1,4,[5],12:e,h:1,5:[R1…])::

    <Subspecies> <O_factor1,O_factor2>:<P1_factor1,P1_factor2>:<P2_factor1,P2_factor2>
    <Subspecies> <O_factor1,O_factor2>:<P1_factor1,P1_factor2>:<P2_factor1,P2_factor2>:<otherH_factor1,otherH_factor2>

- Missing antigens should be specified using '–' (e.g. I 4,12,27:b:– or I 1,9,12:–:–). 

- Optional, exclusive, and weakly agglutinable factors should be designated as follows::

    optional            '[]'
    exclusive           '{}'
    weakly agglutinable '()'
    
- Input may contain the terms Nonmotile, Rough, or Mucoid, which will be converted into the appropriate format::

    Nonmotile  :-:-
    Rough      -:
    Mucoid     -:

- Phage conversion factors denoted by underlining in [1]_ are denoted as optional '[]' in SeroTools, with the exception of exclusive phage conversion factors (e.g. {15} and {15,34}). 


.. [1] Grimont PA, Weill FX. Antigenic Formulae of the Salmonella Serovars. 9th. Paris, France: WHO Collaborating Center for Reference and Research on Salmonella, Institut Pasteur; 2007 <https://www.pasteur.fr/sites/default/files/veng_0.pdf>.
