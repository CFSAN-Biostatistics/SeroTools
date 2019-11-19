========
Usage
========

.. highlight:: none

The White-Kauffmann-Le Minor (WKLM) repository can be queried by submitting an appropriate 
serovar name or antigenic formula in an input file composed of a single query per line:: 

    $ serotools query -i <input_file>
    
or as a command line argument::

    $ serotools query -s 'Paratyphi A'
    
Output::

    Input	     Name	      Formula
    Paratyphi A	 Paratyphi A  I [1],2,12:a:[1,5]

SeroTools provides functionality for comparing serovar predictions by evaluating multiple
states of congruence (exact, congruent, minimally congruent, incongruent). Serovar names 
and/or antigenic formulae may be submitted as a tab-delimited input file composed of two 
columns of serovar predictions::  

    serotools compare -i <input_file>

or as command line arguments::

    serotools compare -1 'Hull' -2 'I 16:b:1,2'

Output::

    Subject	 Name  Formula	   Query	   Name	 Formula	 Result
    Hull	 Hull  I 16:b:1,2  I 16:b:1,2  Hull	 I 16:b:1,2	 exact match
      
SeroTools also provides the ability to predict which serovars are most likely
represented by (minimally congruent with) a serovar name or antigenic formula, which
may be submitted in an input file composed of a single query on each line::

    serotools predict -i <input_file>
    
or as a command line argument::

    serotools predict -s 'I 2,12:–:–'

Output::

    Input       MinCon_Name  MinCon_Formula
    I 2,12:–:–  Kiel         I [1],2,12:g,p:–
    I 2,12:–:–  Koessen      I 2,12:l,v:1,5
    I 2,12:–:–  Nitra        I 2,12:g,m:–
    I 2,12:–:–  Paratyphi A  I [1],2,12:a:[1,5]
      
