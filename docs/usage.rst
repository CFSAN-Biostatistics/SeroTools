========
Usage
========

.. highlight:: none

The White-Kauffmann-Le Minor (WKLM) repository may be queried by submitting an appropriate 
serovar name or antigenic formula as a command line argument or in an input file composed 
of a single query per line. 

    serotools query -i <input_file>

    serotools query -s 'Paratyphi A'
    
    
    Input	     Name	      Formula
    Paratyphi A	 Paratyphi A  I [1],2,12:a:[1,5]

SeroTools provides functionality for comparing serovar predictions by evaluating multiple
states of congruence (exact, congruent, minimally congruent, incongruent). Serovar names 
and/or antigenic formulae may be submitted as as a command line argument or in a 
tab-delimited input file composed of two columns of serovar predictions for comparison.  

    serotools compare -i <input_file>

    serotools compare -1 'Hull' -2 'I 16:b:1,2'

    
    Subject	 Name  Formula	   Query	   Name	 Formula	 Result
    Hull	 Hull  I 16:b:1,2  I 16:b:1,2  Hull	 I 16:b:1,2	 exact match
      
SeroTools also provides the ability to predict which serovars are most likely
represented by (minimally congruent with) a serovar name or antigenic formula, which
may be submitted as a command line argument or in an input file composed of a single 
query on each line.

    serotools predict -i <input_file>

    serotools predict -s 'I 2,12:–:–'

    
    Input       MinCon_Name  MinCon_Formula
    I 2,12:–:–  Kiel         I [1],2,12:g,p:–
    I 2,12:–:–  Koessen      I 2,12:l,v:1,5
    I 2,12:–:–  Nitra        I 2,12:g,m:–
    I 2,12:–:–  Paratyphi A  I [1],2,12:a:[1,5]
      