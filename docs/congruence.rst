===========
Congruence
===========

.. highlight:: none

SeroTools evaluates multiple levels of congruence for comparisons between serovar designations. 


exact
-----

Exact matches must meet one of the following criteria:

- Two serovar designations are the identical string::

    Corvallis                  Corvallis
    I 8,[20]:z4,z23:[z6]       I 8,[20]:z4,z23:[z6]
- Every antigenic factor (**required** or **optional**) matches::

    Corvallis                  I 8,[20]:z4,z23:[z6]
    I 8,[20]:z4,z23:[z6]       I 8,20:z4,z23:z6
    I 1,3,10,19:f,g,t:1,(2),7  I 1,3,10,19:f,g,t:1,2,7
- Neither serovar designation includes any antigenic factors, and the subspecies designations match::

    I ::                       I –:–:–
    II :                       II –:  


congruent
---------

Congruent matches must meet the following criteria:

- The subspecies field must be present for both serovars or neither.

- All **required** antigenic factors match. For example::

    I 6,7,14:g,m,s:–          I 6,7,[14],[54]:g,m,[p],s:–
    I 6,7:g,m,s:–             I 6,7,[14],[54]:g,m,[p],s:[1,2,7]
    Amager var. 15+           Amager
    I 3,15:y:1,2:[z45]        I 3,{10}{15}:y:1,2:[z45]
    6,7:k:[z6]                6,7:k:–                       


minimally congruent
-------------------

Minimally congruent matches must meet the following criteria:

- Every antigen of at least one serovar can be considered a formal subset of the corresponding antigen (no direct conflicts). For example::

    I 6,7,14,[54]:g,m,[p],s:–     6,7,[14],[54]:g,m,[p],s:–
    I                             I 6,7,8,[14],[54]:g,m,[p],s:–
    I 7:g:–                       I 6,7:g,m,s:–
    Gallinarum                    Enteritidis
* Note - the empty set (–) is a subset of every set

The minimally congruent designation is unique to SeroTools and is useful for distinguishing between two scenarios: 

- Serovars which differ due to sample misannotation (incongruent)

- Serovars derived from correctly annotated samples with variation based solely on missing information. When comparing serovar designations, minor differences may be expected due to method-specific irregularities, for example reagent variation for laboratory-based techniques or the presence of nonproductive genomic data when comparing antigenic agglutination to *in silico*-based techniques. Our assumption is that these minor method-specific differences are more likely manifested as missing data (e.g. all but one of the correct factors were detected) than direct conflicts. 


incongruent
-----------
Any comparison which is not minimally congruent. For example::

    I                             II     
    I 1:                          1 2:
    Javiana                       Saintpaul
    I 7,8:g,m,s:–                 I 6,7,[14],[54]:g,m,[p],s:[1,2,7]
    I 4,5:a,b:6,7                 I 5:a,b,c:6,7
