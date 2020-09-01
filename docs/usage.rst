========
Usage
========

.. highlight:: none

SeroTools provides methods for querying and comparing serovar names and antigenic formulae, 
as well as determining the most abundant serovar for a cluster of isolates.

.. _query-label:

query
-----

Query the White-Kauffmann-Le Minor (WKL) repository by submitting one of more 
serovar names or antigenic formulas in an input file composed of a single query per line:: 

    $ serotools query -i <input_file>
    
or as a command line argument::

    $ serotools query -s 'Paratyphi A'
    
Output::

    Input        Name         Formula             Match
    Paratyphi A  Paratyphi A  I [1],2,12:a:[1,5]  exact

.. _compare-label:

compare
-------

Compare serovar predictions by evaluating multiple states of congruence (exact, congruent,
minimally congruent, incongruent). Serovar names and/or antigenic formulae may be submitted 
in a tab-delimited input file composed of two columns of serovar predictions::  

    $ serotools compare -i <input_file>

or as command line arguments::

    $ serotools compare -1 'Hull' -2 'I 16:b:1,2'

Output::

    Serovar1    Name    Formula     Serovar2    Name    Formula     Result
    Hull        Hull    I 16:b:1,2  I 16:b:1,2  Hull    I 16:b:1,2  exact

.. _cluster-label:

cluster
-------
Determine the most abundant serovar(s) for one or more clusters of isolates. Input data must be 
submitted in the form of a tab-delimited file in which each line consists of a cluster ID and one serovar as follows::

Input File - example.txt::

    cluster1	Dunkwa
    cluster1	Dunkwa
    cluster1	Utah
    cluster2	Hull
    
::

    $ serotools cluster -i example.txt
    
Output::

    ClusterID   ClusterSize Input   Name    Formula      P_Exact  P_Congruent	P_MinCon
    cluster1    2           Dunkwa  Dunkwa  I 6,8:d:1,7  0.6667   0.6667        0.6667
    cluster2    1           Hull    Hull    I 16:b:1,2   1.0      1.0           1.0
    
