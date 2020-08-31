===========
Repository
===========

.. highlight:: none

SeroTools provides a repository of the White-Kauffmann-Le Minor (WKL) Salmonella serotyping scheme based on these :ref:`references-label` in the following formats:

- Python data structures (`serotools.py <https://github.com/CFSAN-Biostatistics/SeroTools/blob/master/serotools/serotools.py>`__)

  - pandas DataFrame:: 
  
      wklm_df
    
  - Dictionaries::
  
      wklm_name_to_formula
      wklm_formula_to_name
    
  - Lists with common indexing::
  
      wklm_name
      std_wklm_name (standardized for matching)
      wklm_formula
      std_wklm_formula (standardized for matching)
      wklm_sp (species)
      wklm_subsp (subspecies)
      wklm_O
      wklm_P1
      wklm_P2
      wklm_other_H
      wklm_group (O group)
      wklm_old_group (previous O group)
    
- An Excel spreadsheet (`White-Kauffman-LeMinor-Scheme.xlsx <https://github.com/CFSAN-Biostatistics/SeroTools/blob/master/wklm_scheme/White-Kauffman-LeMinor-Scheme.xlsx>`__)

- A tab-delimited text file (`White-Kauffman-LeMinor_scheme.tsv <https://github.com/CFSAN-Biostatistics/SeroTools/blob/master/wklm_scheme/White-Kauffman-LeMinor_scheme.tsv>`__)
