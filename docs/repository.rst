===========
Repository
===========

.. highlight:: none

SeroTools provides a repository of the White-Kauffmann-Le Minor (WKLM) Salmonella serotyping scheme based on these `references <references.rst>`__ in the following formats:

- Python data structures (`serotools.py <serotools/serotools.py>`__)

  - pandas DataFrame 
  
    - wklm_df
    
  - Dictionaries
  
    - wklm_name_to_formula
    
    - wklm_formula_to_name
    
  - Lists with common indexing
  
    - wklm_name
    
    - std_wklm_name (standardized for matching)
    
    - wklm_formula
    
    - std_wklm_formula (standardized for matching)
    
    - wklm_sp (species)
    
    - wklm_subsp (subspecies)
    
    - wklm_O
    
    - wklm_P1
    
    - wklm_P2
    
    - wklm_other_H
    
    - wklm_group
    
    - wklm_old_group
    
- An Excel spreadsheet (`White-Kauffman-LeMinor-Scheme.xlsx <wklm_scheme/White-Kauffman-LeMinor-Scheme.xlsx>`__)

- A tab-delimited text file (`White-Kauffman-LeMinor_scheme.tsv <wklm_scheme/White-Kauffman-LeMinor_scheme.tsv>`__)
