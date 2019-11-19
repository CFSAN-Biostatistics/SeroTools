#!/usr/bin/env python3

import pytest
import sys
from serotools import serotools as sero
from serotools.serotools import WKLMSerovar, SeroComp, InvalidInput

def test_wklm_repository():

    """WKLM serovar dicts and lists should be the same length""" 
    # sero.sero_formula_to_name is shorter due to names with identical formulas  
    lst_len = len(sero.sero_name)
    assert all(len(lst) == lst_len for lst in [sero.std_sero_name, 
            sero.sero_formula, sero.std_sero_formula, sero.sero_subsp, sero.sero_O, 
            sero.sero_P1, sero.sero_P2, sero.sero_other_H, sero.sero_name_to_formula]) == True   
     

def test_WKLMSerovar():

    input_output = {'Paratyphi A': ['Paratyphi A','Paratyphi A','I [1],2,12:a:[1,5]','I','[1],2,12','a','[1,5]',''],
                    'Westhampton var. 15+,34+': ['Westhampton var. 15+,34+','Westhampton var. 15+,34+','I 3,15,34:g,s,t:–:[z37]','I','3,15,34','g,s,t','–','[z37]'],
                    'westhampton var. 15+,34+': ['westhampton var. 15+,34+','Westhampton var. 15+,34+','I 3,15,34:g,s,t:–:[z37]','I','3,15,34','g,s,t','–','[z37]'],
                    'I 1,4,[5],12,27:z4,z23:[1,2]': ['I 1,4,[5],12,27:z4,z23:[1,2]', 'Stanleyville var. 27+', 'I 1,4,[5],12,27:z4,z23:[1,2]', 'I', '1,4,[5],12,27', 'z4,z23', '[1,2]', ''], 
                    'IV Houten': ['IV Houten','IV 43:z4,z23:–','IV 43:z4,z23:–','IV','43','z4,z23','–',''],
                    'IV houten': ['IV houten','IV 43:z4,z23:–','IV 43:z4,z23:–','IV','43','z4,z23','–',''],
                    'Houten': ['Houten','IV 43:z4,z23:–','IV 43:z4,z23:–','IV','43','z4,z23','–',''],
                    'Illinois': ['Illinois','Lexington var. 15+,34+','I 3,15,34:z10:1,5:[z49]','I','3,15,34','z10','1,5','[z49]'],
                    'Java': ['Java','Paratyphi B var. L(+) tartrate (= d–tartrate)+','I [1],4,[5],12:b:1,2:[z5],[z33]','I','[1],4,[5],12','b','1,2','[z5],[z33]'],
                    'Paratyphi B l(+)': ['Paratyphi B l(+)','Paratyphi B var. L(+) tartrate (= d–tartrate)+','I [1],4,[5],12:b:1,2:[z5],[z33]','I','[1],4,[5],12','b','1,2','[z5],[z33]'],
                    'I [1],2,12:a:[1,5]': ['I [1],2,12:a:[1,5]','Paratyphi A','I [1],2,12:a:[1,5]','I','[1],2,12','a','[1,5]',''],
                    'I 67:r:1,2': ['I 67:r:1,2','Crossness','I 67:r:1,2','I','67','r','1,2',''],
                    'i 67:R:1,2': ['i 67:R:1,2','Crossness','I 67:r:1,2','I','67','r','1,2',''],
                    'IV': ['IV','NA','IV','IV','','','','']}
    
    """No arguments"""
    with pytest.raises(TypeError):
        obj = WKLMSerovar()
                
    """Invalid input - empty"""
    with pytest.raises(InvalidInput):
        empty_obj = WKLMSerovar('')

    """Invalid input -  space"""
    with pytest.raises(InvalidInput):
        space_obj = WKLMSerovar(' ')
        
    """Invalid serovar name"""
    with pytest.raises(InvalidInput):
        bad_name_obj = WKLMSerovar('Test')

    """Serovar name as input"""
    assert WKLMSerovar('Paratyphi A').get_attributes() == input_output['Paratyphi A']
    
    """Last name in list – check for discrepancies in list indices"""
    assert WKLMSerovar('Westhampton var. 15+,34+').get_attributes() \
        == input_output['Westhampton var. 15+,34+']
    
    """Lowercase name"""
    assert WKLMSerovar('westhampton var. 15+,34+').get_attributes() \
        == input_output['westhampton var. 15+,34+']

    """Check special variant"""
    assert WKLMSerovar('I 1,4,[5],12,27:z4,z23:[1,2]').get_attributes() \
        == input_output['I 1,4,[5],12,27:z4,z23:[1,2]']

    """Default check, withdrawn names in sero_old_to_new"""
    assert WKLMSerovar('IV Houten').get_attributes() == input_output['IV Houten']

    """Lowercase, withdrawn names in sero_old_to_new"""
    assert WKLMSerovar('IV houten').get_attributes() == input_output['IV houten']

    """No subsp, withdrawn names in sero_old_to_new"""
    assert WKLMSerovar('Houten').get_attributes() == input_output['Houten']
                
    """Named variant, withdrawn names in sero_old_to_new"""
    assert WKLMSerovar('Illinois').get_attributes() == input_output['Illinois']
        
    """Named variant is Paratyphi B d–tartrate+, withdrawn names in sero_old_to_new"""
    assert WKLMSerovar('Java').get_attributes() == input_output['Java']

    """Input is variation of Paratyphi B var. L(+) tartrate (= d–tartrate)+"""
    assert WKLMSerovar('Paratyphi B l(+)').get_attributes() == input_output['Paratyphi B l(+)']

    """Antigenic formula as input"""
    assert WKLMSerovar('I [1],2,12:a:[1,5]').get_attributes() == input_output['I [1],2,12:a:[1,5]']
        
    """Last formula in list – check for discrepancies in list indices"""
    assert WKLMSerovar('I 67:r:1,2').get_attributes() == input_output['I 67:r:1,2']

    """Check alphacase insensitive"""
    assert WKLMSerovar('i 67:R:1,2').get_attributes() == input_output['i 67:R:1,2']

    """Subsp only"""
    assert WKLMSerovar('IV').get_attributes() == input_output['IV']


def test_invalid_SeroComp():

    """No arguments"""
    with pytest.raises(TypeError):
        obj = SeroComp()
                
    """Invalid input - both empty"""
    assert SeroComp('', '').result == 'invalid input'

    """Invalid input - empty subj"""
    assert SeroComp('', 'I').result == 'invalid input'

    """Invalid input - empty query"""
    assert SeroComp('I', '').result == 'invalid input'

    """Invalid input -  spaces"""
    assert SeroComp(' ', ' ').result == 'invalid input'
        
    """Invalid input -  space for subj"""
    assert SeroComp(' ', 'I').result == 'invalid input'

    """Invalid input -  space for query"""
    assert SeroComp('I', ' ').result == 'invalid input'

    """Invalid input -  en dash only"""
    assert SeroComp('I', '–').result == 'invalid input'
    
    """Invalid input -  semicolon only"""
    assert SeroComp(':', 'I').result == 'invalid input' # subsp required if no antigens

    """Invalid input -  missing antigens only"""
    assert SeroComp('–:–:–', 'I').result == 'invalid input' # subsp required if no antigens

    """Invalid input -  bad name"""
    assert SeroComp('Test', 'I').result == 'invalid input'

def test_exact_SeroComp():
     
    """Exact match -  subsp only"""
    assert SeroComp('I', 'I').result == 'exact match'

    """Exact match -  subsp only lc"""
    assert SeroComp('I', 'i').result == 'exact match'

    """Exact match -  subsp and missing antigens"""
    assert SeroComp('I –:–:–', 'I').result == 'exact match'

    """Exact match -  subsp only"""
    assert SeroComp('I', 'I').result == 'exact match'

    """Optional factor is present [14] - both still match serovar"""
    assert SeroComp('I 6,7,14,[54]:g,m,[p],s:[1,2,7]', 'I 6,7,[14],[54]:g,m,[p],s:[1,2,7]').result \
        == 'exact match'

    """Name (with whitespace) lowercase"""
    assert SeroComp('paratyphi b', 'Paratyphi B').result == 'exact match'
 
    """Paratyphi B variant"""
    assert SeroComp('Paratyphi B', 'Paratyphi B var. L(+) tartrate (= d–tartrate)+').result \
        == 'exact match'

    """Paratyphi B named variant - Java"""
    assert SeroComp('Paratyphi B', 'Java').result == 'exact match'
       
    """Java to Paratyphi B var. L(+) tartrate (= d–tartrate)+"""
    assert SeroComp('Java', 'Paratyphi B var. L(+) tartrate (= d–tartrate)+').result \
        == 'exact match'
        
    """Paratyphi B variant alternate designation"""
    assert SeroComp('Paratyphi B var. dt(+)', 'Paratyphi B var. L(+) tartrate (= d–tartrate)+').result \
        == 'exact match'
  
    """Compare serovar name to formula"""
    assert SeroComp('Montevideo', 'I 6,7,[14],[54]:g,m,[p],s:[1,2,7]').result == 'exact match'
     
    """Compare serovar name to altered formula - optional factor is present"""
    assert SeroComp('Montevideo', 'I 6,7,14,[54]:g,m,[p],s:[1,2,7]').result == 'exact match'


def test_congruent_SeroComp():

    """Missing subsp, matching formula"""
    assert SeroComp('I 6,7,14,[54]:g,m,[p],s:–', '6,7,[14],[54]:g,m,[p],s:–').result \
        == 'congruent'

    """Optional factor present (14) and optional factor not present ([54])"""
    assert SeroComp('I 6,7,14:g,m,[p],s:–', 'I 6,7,[14],[54]:g,m,[p],s:–').result \
        == 'congruent'

    """All optional factors not present"""
    assert SeroComp('I 6,7:g,m,s:–', 'I 6,7,[14],[54]:g,m,[p],s:[1,2,7]').result == 'congruent'

    """All optional factors not present - input reversed"""
    assert SeroComp('I 6,7,[14],[54]:g,m,[p],s:[1,2,7]','I 6,7:g,m,s:–').result == 'congruent'

    """Serovar name to formula - all optional factors not present"""
    assert SeroComp('Montevideo', 'I 6,7:g,m,s:–').result == 'congruent'

    """Serovar name to named variant with exclusive factors"""
    assert SeroComp('Amager', 'Amager var. 15+').result == 'congruent'

    """Serovar name to formula with exclusive factor"""
    assert SeroComp('Amager', 'I 3,10:y:1,2:[z45]').result == 'congruent'

    """Serovar name to formula with exclusive factor"""
    assert SeroComp('Amager', 'I 3,15:y:1,2:[z45]').result == 'congruent'

    """Serovar name to named variant differing only in presence of optional factor"""
    assert SeroComp('Stanleyville','Stanleyville var. 27+').result == 'congruent'

    """Serovar name to formula differing only in presence of optional factor"""
    assert SeroComp('Stanleyville','I 1,4,[5],12,27:z4,z23:[1,2]').result == 'congruent'


def test_min_congruent_SeroComp():

    """Missing or extra O factor (6)"""
    assert SeroComp('I 7:g,m,s:–', 'I 6,7:g,m,s:–').result == 'minimally congruent'
  
    """Missing or extra O factor (6) and missing P1 factors (m,s)"""
    assert SeroComp('I 7:g:–', 'I 6,7:g,m,s:–').result == 'minimally congruent'
   
    """Subsp only to formula"""
    assert SeroComp('I', 'I 6,7,8,[14],[54]:g,m,[p],s:–').result == 'minimally congruent'
   
    """Serovar name compared to formula missing required factor"""
    assert SeroComp('Montevideo', 'I 7:g,m,s:–').result == 'minimally congruent'

    """Serovar name compared to formula missing extra factor"""
    assert SeroComp('Montevideo', 'I 6,7,8,[14],[54]:g,m,[p],s:[1,2,7]').result == 'minimally congruent'

    """Serovar name compared to formula missing extra factor - reverse"""
    assert SeroComp('I 6,7,8,[14],[54]:g,m,[p],s:[1,2,7]', 'Montevideo').result == 'minimally congruent'

    """Gallinarum (I 1,9,12:–:–) to Enteritidis (I 1,9,12:g,m:–)"""
    assert SeroComp('Gallinarum', 'Enteritidis').result == 'minimally congruent'


def test_incongruent_SeroComp():

    assert SeroComp('I', 'II').result == 'incongruent'
    assert SeroComp('I 1:', 'I 2:').result == 'incongruent'  
    assert SeroComp('Javiana', 'Saintpaul').result == 'incongruent'
   
    """Serovar name compared to formula missing required factor and extra factor present"""
    assert SeroComp('Montevideo', 'I 7,8:g,m,s:–').result == 'incongruent'


def test_compare(tmpdir,capsys):                        

    """Input file"""    
    file1 = tmpdir.join('wklm_compare.tsv')
    with open(str(file1), 'w') as f1:
        f1.write('{}\t{}\n'.format('Kumasi','I 30:z10:e,n,z15'))
    
    sero.compare(input_file=file1)
    in_cap = capsys.readouterr()
    in_expected = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                  'Subject','Name','Formula','Query','Name','Formula','Result',
                  'Kumasi','Kumasi','I 30:z10:e,n,z15','I 30:z10:e,n,z15','Kumasi','I 30:z10:e,n,z15','exact match')      
    assert in_cap.out == in_expected

    """Input file - with header"""
    file2 = tmpdir.join('wklm_compare_header.tsv')
    with open(str(file2), 'w') as f2:
        f2.write('{}\t{}\n'.format('Serovar1','Serovar2'))
        f2.write('{}\t{}\n'.format('Kumasi','I 30:z10:e,n,z15'))

    sero.compare(input_file=file2, header=True)
    h_in_cap = capsys.readouterr()
    h_in_expected = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                  'Serovar1','Name','Formula','Serovar2','Name','Formula','Result',
                  'Kumasi','Kumasi','I 30:z10:e,n,z15','I 30:z10:e,n,z15','Kumasi','I 30:z10:e,n,z15','exact match')      
    assert h_in_cap.out == h_in_expected
    
    """Command line args"""
    sero.compare(subj='Kumasi',query='I 30:z10:e,n,z15')
    args_cap = capsys.readouterr()    
    args_expected = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                  'Subject','Name','Formula','Query','Name','Formula','Result',
                  'Kumasi','Kumasi','I 30:z10:e,n,z15','I 30:z10:e,n,z15','Kumasi','I 30:z10:e,n,z15','exact match')
    assert args_cap.out == args_expected
 
    """Invalid input""" 
    sero.compare(subj='test',query='I 30:z10:e,n,z15')
    inv_cap = capsys.readouterr()    
    inv_expected = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                  'Subject','Name','Formula','Query','Name','Formula','Result',
                  'test','NA','NA','I 30:z10:e,n,z15','Kumasi','I 30:z10:e,n,z15','invalid input')
    assert inv_cap.out == inv_expected
  

def test_formula_to_fields():

    """Default formula"""
    formula = 'I 3,{10}{15}{15,34}:z10:1,5:[z49]'
    assert sero.formula_to_fields(formula) == ['I','3,{10}{15}{15,34}','z10','1,5','[z49]']

    """No subsp"""
    formula = '3,{10}{15}{15,34}:z10:1,5:[z49]'
    assert sero.formula_to_fields(formula) == ['','3,{10}{15}{15,34}','z10','1,5','[z49]']
        
    """Subsp only"""
    formula = 'I'
    assert sero.formula_to_fields(formula) == ['I','','','','']
       
    """Invalid formula """
    with pytest.raises(InvalidInput):
        sero.formula_to_fields('Enteritidis')
        
        
def test_is_min_subset():

    assert sero.is_min_subset('5','5') == True
    assert sero.is_min_subset('5','4') == False
    
    """z15 is a proper subset of e,n,z15"""    
    assert sero.is_min_subset('e,n,z15','z15') == True
    
    """Neither is a subset"""
    assert sero.is_min_subset('e,n,z15','z39') == False
    
    """e,n,z15 is a proper subset of e,n,z15,z39"""
    assert sero.is_min_subset('e,n,z15','e,n,z15,z39') == True
    
    """Neither is a subset"""
    assert sero.is_min_subset('e,n,z15','e,n,z39') == False
    
    """Both are proper subsets"""
    assert sero.is_min_subset('6,7,8,[14],[54]','6,7,8') == True
    
    """6,7,8 is a proper subset of 6,7,8,9"""
    assert sero.is_min_subset('6,7,8,[14],[54]','6,7,8,9') == True
    
    """Neither is a subset"""    
    assert sero.is_min_subset('6,7,8,10,[14],[54]','6,7,8,9') == False
    
    """Both are proper subsets"""
    assert sero.is_min_subset('g,m,[p],s','g,m,p,s') == True
    
    """m,p,s is a proper subset of g,m,[p],s"""
    assert sero.is_min_subset('g,m,[p],s','m,p,s') == True
    
    """m is a proper subset of g,[m],[s],[t]"""
    assert sero.is_min_subset('g,[m],[s],[t]','m') == True
    
    """the empty set is a proper subset of the empty set"""
    assert sero.is_min_subset('[1,2,7]','–') == True
    
    """the empty set is a proper subset of 1,2,7"""
    assert sero.is_min_subset('1,2,7','–') == True

    """the empty set is a proper subset of the empty set"""    
    assert sero.is_min_subset('–','–') == True

    """the empty set is a proper subset of the empty set"""    
    assert sero.is_min_subset('[1,2,7]','[5]') == True
    
    
def test_is_name():

    assert sero.is_name('Enteritidis') == True
    assert sero.is_name('Paratyphi B') == True
    assert sero.is_name('Paratyphi B var. L(+) tartrate (= d–tartrate)+') == True
    assert sero.is_name('Senftenberg or Dessau') == True
    assert sero.is_name('Senftenberg/Dessau') == True
    assert sero.is_name('IV Houten') == True
    assert sero.is_name('I 1,2,12:a:[1,5]') == False
    assert sero.is_name('I') == False
    assert sero.is_name('I :') == False
    assert sero.is_name('Paratyphi A (I 1,2,12:a:[1,5])') == False
    assert sero.is_name('') == True   # no checks here
    assert sero.is_name(' ') == True  # no checks here


def test_is_opt_factor():

    assert sero.is_opt_factor('g','g,[m],(s),{t}') == False
    assert sero.is_opt_factor('m','g,[m],(s),{t}') == True
    assert sero.is_opt_factor('s','g,[m],(s),{t}') == True
    assert sero.is_opt_factor('t','g,[m],(s),{t}') == True
    assert sero.is_opt_factor('x','g,[m],(s),{t}') == False


def test_is_opt_subset():

    """Any differences must be due to optional factors"""    

    assert sero.is_opt_subset('5','5') == True
    assert sero.is_opt_subset('5','4') == False  
    assert sero.is_opt_subset('e,n,z15','z15') == False
    assert sero.is_opt_subset('e,n,z15','z39') == False
    assert sero.is_opt_subset('e,n,z15','e,n,z15,z39') == False
    assert sero.is_opt_subset('e,n,z15','e,n,z39') == False
    assert sero.is_opt_subset('6,7,8,[14],[54]','6,7,8') == True
    assert sero.is_opt_subset('6,7,8,[14],[54]','6,7,8,9') == False
    assert sero.is_opt_subset('6,7,8,10,[14],[54]','6,7,8,9') == False
    assert sero.is_opt_subset('g,m,[p],s','g,m,p,s') == True
    assert sero.is_opt_subset('g,m,[p],s','g,m,s') == True
    assert sero.is_opt_subset('g,m,[p],s','m,p,s') == False
    assert sero.is_opt_subset('g,[m],[s],[t]','m') == False
    assert sero.is_opt_subset('[1,2,7]','–') == True
    assert sero.is_opt_subset('1,2,7','–') == False
    assert sero.is_opt_subset('–','–') == True
    assert sero.is_opt_subset('[1,2,7]','[5]') == True
'g,m,[p],s,[t]','f,g,m,p,s'

def test_matching_indices():

    assert sero.matching_indices('5',['1,2','2,5,6','10,15','2,5,7,10']) == [1,3]
    assert sero.matching_indices('–',['1,2','2,5,6,z15','–','–']) == [2,3]
    assert sero.matching_indices('i',['I','I','V','III']) == [0,1]

  
def test_max_factors():

    assert sero.max_factors('[e,n,x]') == {'e', 'n', 'x'}
    assert sero.max_factors('[e],(n),{x}') == {'e', 'n', 'x'}
    assert sero.max_factors('1,2,5') == {'1', '2', '5'}
    assert sero.max_factors('[1],2,5') == {'1', '2', '5'}
    assert sero.max_factors('-') == {'-'}

    
def test_min_congruent_serovars(): 

    assert sero.min_congruent_serovars('I 1,9,12:b:–') \
        == {'Onarimon': 'I [1],9,12:b:1,2', 'Frintrop': 'I [1],9,12:b:1,5'}
    assert sero.min_congruent_serovars('I 1,4,[5],12:b:1,2:[z5],[z33]') \
        == {'Paratyphi B var. L(+) tartrate (= d–tartrate)+': 'I [1],4,[5],12:b:1,2:[z5],[z33]',
            'Paratyphi B': 'I [1],4,[5],12:b:1,2:[z5],[z33]'}
    

def test_min_factors():

    assert sero.min_factors('[e,n,x]') == set()
    assert sero.min_factors('[e],(n),{x}') == set()
    assert sero.min_factors('1,2,5') == {'1', '2', '5'}
    assert sero.min_factors('[1],2,5') == {'2', '5'}
    assert sero.min_factors('–') == {'–'}


def test_predict(tmpdir,capsys):                        

    """Input file""" 
    file = tmpdir.join('wklm_predict.tsv')   
    with open(str(file), 'w') as f:
        f.write('{}\n'.format('I 2,12:–:–'))
    
    sero.predict(input_file=file)
    in_cap = capsys.readouterr()
    in_expected = '{}\t{}\t{}\n{}\t{}\t{}\n{}\t{}\t{}\n{}\t{}\t{}\n{}\t{}\t{}\n'.format(
                  'Input','MinCon_Name','MinCon_Formula',
                  'I 2,12:–:–','Kiel','I [1],2,12:g,p:–',
                  'I 2,12:–:–','Koessen','I 2,12:l,v:1,5',
                  'I 2,12:–:–','Nitra','I 2,12:g,m:–',
                  'I 2,12:–:–','Paratyphi A','I [1],2,12:a:[1,5]')
    assert in_cap.out == in_expected
 
    """Input file - output transformed""" 
    sero.predict(input_file=file, transform=True)
    t_cap = capsys.readouterr()
    t_expected = '{}\t{}\n'.format(
                  'I 2,12:–:–', 'Kiel=I [1],2,12:g,p:–;Koessen=I 2,12:l,v:1,5;Nitra=I 2,12:g,m:–;Paratyphi A=I [1],2,12:a:[1,5]')      
    assert t_cap.out == t_expected
   
    """Command line arg"""
    sero.predict(serovar='I 2,12:–:–')
    arg_cap = capsys.readouterr()    
    arg_expected = '{}\t{}\t{}\n{}\t{}\t{}\n{}\t{}\t{}\n{}\t{}\t{}\n{}\t{}\t{}\n'.format(
                  'Input','MinCon_Name','MinCon_Formula',
                  'I 2,12:–:–','Kiel','I [1],2,12:g,p:–',
                  'I 2,12:–:–','Koessen','I 2,12:l,v:1,5',
                  'I 2,12:–:–','Nitra','I 2,12:g,m:–',
                  'I 2,12:–:–','Paratyphi A','I [1],2,12:a:[1,5]')
    assert arg_cap.out == arg_expected

    """Invalid input""" 
    sero.predict(serovar='test')
    inv_cap = capsys.readouterr()    
    inv_expected ='{}\t{}\t{}\n{}\t{}\t{}\n'.format(
                      'Input','MinCon_Name','MinCon_Formula',
                      'test','NA','NA')
    assert inv_cap.out == inv_expected

    """Invalid input - output transformed""" 
    sero.predict(serovar='test',transform=True)
    t_inv_cap = capsys.readouterr()    
    t_inv_expected ='{}\t{}\n'.format(
                        'test','NA')
    assert t_inv_cap.out == t_inv_expected
 
def test_prep():

    assert sero.prep('I 1,(4),[5],{10}{15}{15,34}:b:1,2:[z5],[z33]') == 'i 1,4,5,10,15,15,34:b:1,2:z5,z33'
    assert sero.prep('I [1],4,[5],12:e,h:1,5:[R1…]') == 'i 1,4,5,12:e,h:1,5:r1'
    assert sero.prep('I [1],4,[5],12:e,h:1,5:[R1...]') == 'i 1,4,5,12:e,h:1,5:r1'

    
def test_query(tmpdir,capsys):                        

    file = tmpdir.join('wklm_query.tsv')   
    with open(str(file), 'w') as f:
        f.write('{}\n'.format('Kumasi'))
    
    """Input file""" 
    sero.query(input_file=file)
    in_cap = capsys.readouterr()
    in_expected ='{}\t{}\t{}\n{}\t{}\t{}\n'.format(
                 'Input','Name','Formula',
                 'Kumasi','Kumasi','I 30:z10:e,n,z15')       
    assert in_cap.out == in_expected
    
    """Command line arg"""
    sero.query(serovar='Kumasi')
    arg_cap = capsys.readouterr()    
    arg_expected ='{}\t{}\t{}\n{}\t{}\t{}\n'.format(
                 'Input','Name','Formula',
                 'Kumasi','Kumasi','I 30:z10:e,n,z15')
    assert arg_cap.out == arg_expected

    """Invalid input""" 
    sero.query(serovar='test')
    inv_cap = capsys.readouterr()    
    inv_expected ='{}\t{}\t{}\n{}\t{}\t{}\n'.format(
                 'Input','Name','Formula',
                 'test','NA','NA')
    assert inv_cap.out == inv_expected

def test_split_input():

    """Multiple serovar predictions as input"""
    assert sero.split_input('Paratyphi C or Choleraesuis or Typhisuis*') \
        == ['Paratyphi C', 'Choleraesuis', 'Typhisuis*']
    assert sero.split_input('Tananarive/Brunei') == ['Tananarive', 'Brunei']


def test_standardize_formula():

    """Missing antigens"""
    assert sero.standardize_formula('I 4,5,12:nonmotile') == 'I 4,5,12:–:–'
    assert sero.standardize_formula('I 4,5,12:Non-Motile') == 'I 4,5,12:–:–'
    assert sero.standardize_formula('I Rough:i:1,2') == 'I –:i:1,2'
    assert sero.standardize_formula('I Rough:non-motile') == 'I –:–:–'
    assert sero.standardize_formula('I mucoid:i:1,2') == 'I –:i:1,2'
    assert sero.standardize_formula('I 4,5,12::') == 'I 4,5,12:–:–'
    assert sero.standardize_formula('I 4,5,12:i:undetermined') == 'I 4,5,12:i:–'
    assert sero.standardize_formula('I Undetermined:i:1,2') == 'I –:i:1,2'
    assert sero.standardize_formula('I :i:1,2') == 'I –:i:1,2'
    assert sero.standardize_formula(':i:1,2') == '–:i:1,2'
    
    """Separate exclusive factors"""
    assert sero.standardize_formula('I 1,4,5,{10}{15}{15,34}:b:1,2') == 'I 1,4,5,{10},{15},{15,34}:b:1,2'
    
    """Invalid input"""
    with pytest.raises(InvalidInput):
        sero.standardize_formula('Rough:non-motile')
    with pytest.raises(InvalidInput):
        sero.standardize_formula('–:–:–')
    with pytest.raises(InvalidInput):
        sero.standardize_formula('–:–:–:–')
    

def test_standardize_input():

    assert sero.standardize_input('Poona*') == 'Poona'

    """Remove content in parentheses, except for weakly agglutinable factors."""
    assert sero.standardize_input('Cerro var. 14+ (Siegburg)') == 'Cerro var. 14+'
    assert sero.standardize_input('I 1,3,10,19:f,g,t:1,(2),7') == 'I 1,3,10,19:f,g,t:1,(2),7'
    assert sero.standardize_input('IIIb 16:(k):e,n,x,z15') == 'IIIb 16:(k):e,n,x,z15'
    
    """Transform subsp."""
    assert sero.standardize_input('Enterica 16:i:z6') == 'I 16:i:z6'
    assert sero.standardize_input('salamae 16:d:1,5') == 'II 16:d:1,5'
    assert sero.standardize_input('Arizonae 13,23:z4,z24:–') == 'IIIa 13,23:z4,z24:–'
    assert sero.standardize_input('Diarizonae (6),14:r:z') == 'IIIb (6),14:r:z'
    assert sero.standardize_input('houtenae 6,14:z4,z23:–') == 'IV 6,14:z4,z23:–'
    assert sero.standardize_input('bongori 1,40:z35:–') == 'V 1,40:z35:–'
    assert sero.standardize_input('indica 41:b:1,7') == 'VI 41:b:1,7'
    
    """Remove leading/trailing whitespace"""
    assert sero.standardize_input(' Enterica 16:i:z6 ') == 'I 16:i:z6'
    
 
def test_standardize_unicode():
 
    assert sero.standardize_unicode('I [1],4,[5],12:e,h:1,5:[R1...]') == 'I [1],4,[5],12:e,h:1,5:[R1…]'
    assert sero.standardize_unicode('I 4,5,12:i:-') == 'I 4,5,12:i:–'
