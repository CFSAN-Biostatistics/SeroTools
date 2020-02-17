#!/usr/bin/env python3

import sys
import pytest
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from serotools import serotools as st
from serotools.serotools import WKLMSerovar, SeroClust, SeroComp, InvalidInput


def test_wklm_repository():

    """WKLM serovar dicts and lists should be the same length""" 
    # st.wklm_formula_to_name is shorter due to names with identical formulas  
    lst_len = len(st.wklm_name)
    assert all(len(lst) == lst_len for lst in [st.std_wklm_name, st.wklm_formula, 
            st.std_wklm_formula, st.wklm_subsp, st.wklm_O, st.wklm_P1, st.wklm_P2, 
            st.wklm_other_H, st.wklm_group, st.wklm_old_group, st.wklm_name_to_formula]) == True   
     

def test_WKLMSerovar(caplog):    
    
    sero_index=['Name','Std_Name','Formula','Std_Formula','Species','Subspecies',
            'O','P1','P2','other_H','Group','Old_Group','Input']

    """No arguments"""
    with pytest.raises(TypeError):
        obj = WKLMSerovar()
                
    """Invalid input - empty"""
    caplog.clear()
    assert WKLMSerovar('').meta.equals(pd.Series({ 
          'Name':np.nan,'Std_Name':np.nan,'Formula':np.nan,'Std_Formula':np.nan,
          'Species':np.nan,'Subspecies':np.nan,'O':np.nan,'P1':np.nan,'P2':np.nan,
          'other_H':np.nan,'Group':np.nan,'Old_Group':np.nan,'Input':''},
          index=sero_index))
    assert caplog.records[0].message == "The serovar '' was not recognized."
    
    """Invalid input -  space"""
    caplog.clear()
    assert WKLMSerovar(' ').meta.equals(pd.Series({ 
          'Name':np.nan,'Std_Name':np.nan,'Formula':np.nan,'Std_Formula':np.nan,
          'Species':np.nan,'Subspecies':np.nan,'O':np.nan,'P1':np.nan,'P2':np.nan,
          'other_H':np.nan,'Group':np.nan,'Old_Group':np.nan,'Input':' '},
          index=sero_index))
    assert caplog.records[0].message == "The serovar ' ' was not recognized."
    
    """Invalid serovar name"""
    caplog.clear()
    assert WKLMSerovar('Test').meta.equals(pd.Series({ 
          'Name':np.nan,'Std_Name':np.nan,'Formula':np.nan,'Std_Formula':np.nan,
          'Species':np.nan,'Subspecies':np.nan,'O':np.nan,'P1':np.nan,'P2':np.nan,
          'other_H':np.nan,'Group':np.nan,'Old_Group':np.nan,'Input':'Test'},
          index=sero_index))
    assert caplog.records[0].message == "The serovar 'Test' was not recognized."
        
    """Check input attribute"""
    assert WKLMSerovar('Paratyphi A').input == 'Paratyphi A'

    """Check name attribute"""
    assert WKLMSerovar('Paratyphi A').name == 'Paratyphi A'

    """Check formula attribute"""
    assert WKLMSerovar('Paratyphi A').formula == 'I [1],2,12:a:[1,5]'
    
    """Check meta attribute"""
    assert WKLMSerovar('Paratyphi A').meta.equals(pd.Series({ 
          'Name':'Paratyphi A','Std_Name':'paratyphi a',
          'Formula':'I [1],2,12:a:[1,5]','Std_Formula':'i 1,2,12:a:1,5',
          'Species':'enterica','Subspecies':'I','O':'[1],2,12','P1':'a','P2':'[1,5]',
          'other_H':'','Group':'O:2','Old_Group':'A','Input':'Paratyphi A'},
          index=sero_index))
    
    """Last name in list – check for discrepancies in list indices"""
    assert WKLMSerovar('Westhampton var. 15+,34+').meta.equals(pd.Series({ 
        'Name':'Westhampton var. 15+,34+','Std_Name':'westhampton var. 15+,34+',
        'Formula':'I 3,15,34:g,s,t:–:[z37]','Std_Formula':'i 3,15,34:g,s,t:–:z37',
        'Species':'enterica','Subspecies':'I','O':'3,15,34','P1':'g,s,t','P2':'–',
        'other_H':'[z37]','Group':'O:3,10','Old_Group':'E1','Input':'Westhampton var. 15+,34+'},
          index=sero_index))
        
    """Lowercase name"""
    assert WKLMSerovar('westhampton var. 15+,34+').meta.equals(pd.Series({ 
        'Name':'Westhampton var. 15+,34+','Std_Name':'westhampton var. 15+,34+',
        'Formula':'I 3,15,34:g,s,t:–:[z37]','Std_Formula':'i 3,15,34:g,s,t:–:z37',
        'Species':'enterica','Subspecies':'I','O':'3,15,34','P1':'g,s,t','P2':'–',
        'other_H':'[z37]','Group':'O:3,10','Old_Group':'E1','Input':'westhampton var. 15+,34+'},
          index=sero_index))

    """Check special variant"""
    assert WKLMSerovar('I 1,4,[5],12,27:z4,z23:[1,2]').meta.equals(pd.Series({ 
        'Name':'Stanleyville var. 27+','Std_Name':'stanleyville var. 27+',
        'Formula':'I 1,4,[5],12,27:z4,z23:[1,2]','Std_Formula':'i 1,4,5,12,27:z4,z23:1,2',
        'Species':'enterica','Subspecies':'I','O':'1,4,[5],12,27','P1':'z4,z23','P2':'[1,2]',
        'other_H':'','Group':'O:4','Old_Group':'B','Input':'I 1,4,[5],12,27:z4,z23:[1,2]'},
          index=sero_index))

    """Default check, withdrawn names in wklm_old_to_new"""
    assert WKLMSerovar('IV Houten').meta.equals(pd.Series({ 
        'Name':'IV 43:z4,z23:–','Std_Name':'iv 43:z4,z23:–',
        'Formula':'IV 43:z4,z23:–','Std_Formula':'iv 43:z4,z23:–',
        'Species':'enterica','Subspecies':'IV','O':'43','P1':'z4,z23','P2':'–',
        'other_H':'','Group':'O:43','Old_Group':'U','Input':'IV Houten'},
          index=sero_index))

    """Lowercase, withdrawn names in wklm_old_to_new"""
    assert WKLMSerovar('IV houten').meta.equals(pd.Series({ 
        'Name':'IV 43:z4,z23:–','Std_Name':'iv 43:z4,z23:–',
        'Formula':'IV 43:z4,z23:–','Std_Formula':'iv 43:z4,z23:–',
        'Species':'enterica','Subspecies':'IV','O':'43','P1':'z4,z23','P2':'–',
        'other_H':'','Group':'O:43','Old_Group':'U','Input':'IV houten'},
          index=sero_index))

    """No subsp, withdrawn names in wklm_old_to_new"""
    assert WKLMSerovar('Houten').meta.equals(pd.Series({ 
        'Name':'IV 43:z4,z23:–','Std_Name':'iv 43:z4,z23:–',
        'Formula':'IV 43:z4,z23:–','Std_Formula':'iv 43:z4,z23:–',
        'Species':'enterica','Subspecies':'IV','O':'43','P1':'z4,z23','P2':'–',
        'other_H':'','Group':'O:43','Old_Group':'U','Input':'Houten'},
          index=sero_index))
                
    """Named variant, withdrawn names in wklm_old_to_new"""
    assert WKLMSerovar('Illinois').meta.equals(pd.Series({ 
        'Name':'Lexington var. 15+,34+','Std_Name':'lexington var. 15+,34+',
        'Formula':'I 3,15,34:z10:1,5:[z49]','Std_Formula':'i 3,15,34:z10:1,5:z49',
        'Species':'enterica','Subspecies':'I','O':'3,15,34','P1':'z10','P2':'1,5',
        'other_H':'[z49]','Group':'O:3,10','Old_Group':'E1','Input':'Illinois'},
          index=sero_index))
                
    """Named variant is Paratyphi B d–tartrate+, withdrawn names in wklm_old_to_new"""
    assert WKLMSerovar('Java').meta.equals(pd.Series({ 
        'Name':'Paratyphi B var. L(+) tartrate (= d–tartrate)+',
        'Std_Name':'paratyphi b var. l+ tartrate = d–tartrate+',
        'Formula':'I [1],4,[5],12:b:1,2:[z5],[z33]','Std_Formula':'i 1,4,5,12:b:1,2:z5,z33',
        'Species':'enterica','Subspecies':'I','O':'[1],4,[5],12','P1':'b','P2':'1,2',
        'other_H':'[z5],[z33]','Group':'O:4','Old_Group':'B','Input':'Java'},
          index=sero_index))

    """Input is variation of Paratyphi B var. L(+) tartrate (= d–tartrate)+"""
    assert WKLMSerovar('Paratyphi B l(+)').meta.equals(pd.Series({ 
        'Name':'Paratyphi B var. L(+) tartrate (= d–tartrate)+',
        'Std_Name':'paratyphi b var. l+ tartrate = d–tartrate+',
        'Formula':'I [1],4,[5],12:b:1,2:[z5],[z33]','Std_Formula':'i 1,4,5,12:b:1,2:z5,z33',
        'Species':'enterica','Subspecies':'I','O':'[1],4,[5],12','P1':'b','P2':'1,2',
        'other_H':'[z5],[z33]','Group':'O:4','Old_Group':'B','Input':'Paratyphi B l(+)'},
          index=sero_index))

    """Antigenic formula as input"""
    assert WKLMSerovar('I [1],2,12:a:[1,5]').meta.equals(pd.Series({ 
          'Name':'Paratyphi A','Std_Name':'paratyphi a',
          'Formula':'I [1],2,12:a:[1,5]','Std_Formula':'i 1,2,12:a:1,5',
          'Species':'enterica','Subspecies':'I','O':'[1],2,12','P1':'a','P2':'[1,5]',
          'other_H':'','Group':'O:2','Old_Group':'A','Input':'I [1],2,12:a:[1,5]'},
          index=sero_index))
        
    """Last formula in list – check for discrepancies in list indices"""
    assert WKLMSerovar('I 67:r:1,2').meta.equals(pd.Series({ 
          'Name':'Crossness','Std_Name':'crossness',
          'Formula':'I 67:r:1,2','Std_Formula':'i 67:r:1,2',
          'Species':'enterica','Subspecies':'I','O':'67','P1':'r','P2':'1,2',
          'other_H':'','Group':'O:67','Old_Group':'–','Input':'I 67:r:1,2'},
          index=sero_index))

    """Check alphacase insensitive"""
    assert WKLMSerovar('i 67:R:1,2').meta.equals(pd.Series({ 
          'Name':'Crossness','Std_Name':'crossness',
          'Formula':'I 67:r:1,2','Std_Formula':'i 67:r:1,2',
          'Species':'enterica','Subspecies':'I','O':'67','P1':'r','P2':'1,2',
          'other_H':'','Group':'O:67','Old_Group':'–','Input':'i 67:R:1,2'},
          index=sero_index))

    """Subsp only"""
    assert WKLMSerovar('IV').meta.equals(pd.Series({ 
          'Name':np.nan,'Std_Name':np.nan,
          'Formula':'IV','Std_Formula':'iv',
          'Species':'enterica','Subspecies':'IV','O':'–','P1':'–','P2':'–',
          'other_H':'','Group':np.nan,'Old_Group':np.nan,'Input':'IV'},
          index=sero_index))

    """Identical formulas - by name (Choleraesuis and Typhisuis)"""
    assert WKLMSerovar('Choleraesuis').meta.equals(pd.Series({ 
          'Name':'Choleraesuis','Std_Name':'choleraesuis',
          'Formula':'I 6,7:c:1,5','Std_Formula':'i 6,7:c:1,5',
          'Species':'enterica','Subspecies':'I','O':'6,7','P1':'c','P2':'1,5',
          'other_H':'','Group':'O:7','Old_Group':'C1','Input':'Choleraesuis'},
          index=sero_index))

    """Identical formulas - by name (Choleraesuis and Typhisuis)"""
    assert WKLMSerovar('Typhisuis').meta.equals(pd.Series({ 
          'Name':'Typhisuis','Std_Name':'typhisuis',
          'Formula':'I 6,7:c:1,5','Std_Formula':'i 6,7:c:1,5',
          'Species':'enterica','Subspecies':'I','O':'6,7','P1':'c','P2':'1,5',
          'other_H':'','Group':'O:7','Old_Group':'C1','Input':'Typhisuis'},
          index=sero_index))

    """Identical formulas - formula"""
    assert WKLMSerovar('I 6,7:c:1,5').meta.equals(pd.Series({ 
          'Name':'Choleraesuis or Typhisuis','Std_Name':'choleraesuis or typhisuis',
          'Formula':'I 6,7:c:1,5','Std_Formula':'i 6,7:c:1,5',
          'Species':'enterica','Subspecies':'I','O':'6,7','P1':'c','P2':'1,5',
          'other_H':'','Group':'O:7','Old_Group':'C1','Input':'I 6,7:c:1,5'},
          index=sero_index))


def test_invalid_SeroComp():

    """No arguments"""
    with pytest.raises(TypeError):
        obj = SeroComp()
                
    """Invalid input - both empty"""
    assert SeroComp(WKLMSerovar(''),WKLMSerovar('')).result == 'invalid input'

    """Invalid input - empty subj"""
    assert SeroComp(WKLMSerovar(''),WKLMSerovar('I')).result == 'invalid input'

    """Invalid input - empty query"""
    assert SeroComp(WKLMSerovar('I'),WKLMSerovar('')).result == 'invalid input'

    """Invalid input -  spaces"""
    assert SeroComp(WKLMSerovar(' '),WKLMSerovar(' ')).result == 'invalid input'
        
    """Invalid input -  space for subj"""
    assert SeroComp(WKLMSerovar(' '),WKLMSerovar('I')).result == 'invalid input'

    """Invalid input -  space for query"""
    assert SeroComp(WKLMSerovar('I'),WKLMSerovar(' ')).result == 'invalid input'

    """Invalid input -  en dash only"""
    assert SeroComp(WKLMSerovar('I'),WKLMSerovar('–')).result == 'invalid input'
    
    """Invalid input -  semicolon only"""
    assert SeroComp(WKLMSerovar(':'),WKLMSerovar('I')).result == 'invalid input' # subsp required if no antigens

    """Invalid input -  missing antigens only"""
    assert SeroComp(WKLMSerovar('–:–:–'),WKLMSerovar('I')).result == 'invalid input' # subsp required if no antigens

    """Invalid input -  bad name"""
    assert SeroComp(WKLMSerovar('Test'),WKLMSerovar('I')).result == 'invalid input'


def test_exact_SeroComp():
     
    """Exact match -  subsp only"""
    assert SeroComp(WKLMSerovar('I'),WKLMSerovar('I')).result == 'exact'

    """Exact match -  subsp only lc"""
    assert SeroComp(WKLMSerovar('I'),WKLMSerovar('i')).result == 'exact'

    """Exact match -  subsp and missing antigens"""
    assert SeroComp(WKLMSerovar('I –:–:–'),WKLMSerovar('I')).result == 'exact'

    """Exact match -  subsp only"""
    assert SeroComp(WKLMSerovar('I'),WKLMSerovar('I')).result == 'exact'

    """Optional factor is present [14] - both still match serovar"""
    assert SeroComp(WKLMSerovar('I 6,7,14,[54]:g,m,[p],s:[1,2,7]'),
                    WKLMSerovar('I 6,7,[14],[54]:g,m,[p],s:[1,2,7]')).result == 'exact'

    """Name (with whitespace) lowercase"""
    assert SeroComp(WKLMSerovar('paratyphi b'),WKLMSerovar('Paratyphi B')).result == 'exact'
 
    """Paratyphi B variant"""
    assert SeroComp(WKLMSerovar('Paratyphi B'),
                    WKLMSerovar('Paratyphi B var. L(+) tartrate (= d–tartrate)+')).result \
                    == 'exact'

    """Paratyphi B named variant - Java"""
    assert SeroComp(WKLMSerovar('Paratyphi B'),WKLMSerovar('Java')).result == 'exact'
       
    """Java to Paratyphi B var. L(+) tartrate (= d–tartrate)+"""
    assert SeroComp(WKLMSerovar('Java'),
                    WKLMSerovar('Paratyphi B var. L(+) tartrate (= d–tartrate)+')).result \
                    == 'exact'
        
    """Paratyphi B variant alternate designation"""
    assert SeroComp(WKLMSerovar('Paratyphi B var. dt(+)'),
                    WKLMSerovar('Paratyphi B var. L(+) tartrate (= d–tartrate)+')).result \
                    == 'exact'
  
    """Compare serovar name to formula"""
    assert SeroComp(WKLMSerovar('Montevideo'),
                    WKLMSerovar('I 6,7,[14],[54]:g,m,[p],s:[1,2,7]')).result == 'exact'
     
    """Compare serovar name to altered formula - optional factor is present"""
    assert SeroComp(WKLMSerovar('Montevideo'),
                    WKLMSerovar('I 6,7,14,[54]:g,m,[p],s:[1,2,7]')).result == 'exact'


def test_congruent_SeroComp():

    """Optional factor present (14) and optional factor not present ([54])"""
    assert SeroComp(WKLMSerovar('I 6,7,14:g,m,[p],s:–'),
                    WKLMSerovar('I 6,7,[14],[54]:g,m,[p],s:–')).result == 'congruent'

    """All optional factors not present"""
    assert SeroComp(WKLMSerovar('I 6,7:g,m,s:–'),
                    WKLMSerovar('I 6,7,[14],[54]:g,m,[p],s:[1,2,7]')).result == 'congruent'

    """All optional factors not present - input reversed"""
    assert SeroComp(WKLMSerovar('I 6,7,[14],[54]:g,m,[p],s:[1,2,7]'),
                    WKLMSerovar('I 6,7:g,m,s:–')).result == 'congruent'

    """Serovar name to formula - all optional factors not present"""
    assert SeroComp(WKLMSerovar('Montevideo'),
                    WKLMSerovar('I 6,7:g,m,s:–')).result == 'congruent'

    """Serovar name to named variant with exclusive factors"""
    assert SeroComp(WKLMSerovar('Amager'),
                    WKLMSerovar('Amager var. 15+')).result == 'congruent'

    """Serovar name to formula with exclusive factor"""
    assert SeroComp(WKLMSerovar('Amager'),
                    WKLMSerovar('I 3,10:y:1,2:[z45]')).result == 'congruent'

    """Serovar name to formula with exclusive factor"""
    assert SeroComp(WKLMSerovar('Amager'),
                    WKLMSerovar('I 3,15:y:1,2:[z45]')).result == 'congruent'

    """Serovar name to named variant differing only in presence of optional factor"""
    assert SeroComp(WKLMSerovar('Stanleyville'),
                    WKLMSerovar('Stanleyville var. 27+')).result == 'congruent'

    """Serovar name to formula differing only in presence of optional factor"""
    assert SeroComp(WKLMSerovar('Stanleyville'),
                    WKLMSerovar('I 1,4,[5],12,27:z4,z23:[1,2]')).result == 'congruent'


def test_min_congruent_SeroComp():

    """Missing subsp, matching formula"""
    assert SeroComp(WKLMSerovar('I 6,7,14,[54]:g,m,[p],s:–'),
                    WKLMSerovar('6,7,[14],[54]:g,m,[p],s:–')).result == 'minimally congruent'
                    
    """Missing or extra O factor (6)"""
    assert SeroComp(WKLMSerovar('I 7:g,m,s:–'),
                    WKLMSerovar('I 6,7:g,m,s:–')).result == 'minimally congruent'
  
    """Missing or extra O factor (6) and missing P1 factors (m,s)"""
    assert SeroComp(WKLMSerovar('I 7:g:–'),
                    WKLMSerovar('I 6,7:g,m,s:–')).result == 'minimally congruent'
   
    """Subsp only to formula"""
    assert SeroComp(WKLMSerovar('I'),
                    WKLMSerovar('I 6,7,8,[14],[54]:g,m,[p],s:–')).result == 'minimally congruent'
   
    """Serovar name compared to formula missing required factor"""
    assert SeroComp(WKLMSerovar('Montevideo'),
                    WKLMSerovar('I 7:g,m,s:–')).result == 'minimally congruent'

    """Serovar name compared to formula missing extra factor"""
    assert SeroComp(WKLMSerovar('Montevideo'),
                    WKLMSerovar('I 6,7,8,[14],[54]:g,m,[p],s:[1,2,7]')).result == 'minimally congruent'

    """Serovar name compared to formula missing extra factor - reverse"""
    assert SeroComp(WKLMSerovar('I 6,7,8,[14],[54]:g,m,[p],s:[1,2,7]'),
                    WKLMSerovar('Montevideo')).result == 'minimally congruent'

    """Gallinarum (I 1,9,12:–:–) to Enteritidis (I 1,9,12:g,m:–)"""
    assert SeroComp(WKLMSerovar('Gallinarum'),
                    WKLMSerovar('Enteritidis')).result == 'minimally congruent'


def test_incongruent_SeroComp():

    assert SeroComp(WKLMSerovar('I'),WKLMSerovar('II')).result == 'incongruent'
    assert SeroComp(WKLMSerovar('I 1:'),WKLMSerovar('I 2:')).result == 'incongruent'  
    assert SeroComp(WKLMSerovar('Javiana'),WKLMSerovar('Saintpaul')).result == 'incongruent'
   
    """Serovar name compared to formula missing required factor and extra factor present"""
    assert SeroComp(WKLMSerovar('Montevideo'),WKLMSerovar('I 7,8:g,m,s:–')).result == 'incongruent'

    """Individual factors are subsets, but neither formula is a proper subset of the other"""
    assert SeroComp(WKLMSerovar('I 4,5:a,b:6,7'),WKLMSerovar('I 5:a,b,c:6,7')).result == 'incongruent'


def test_SeroClust():
   
    results_cols=['ClusterID','ClusterSize','Input','Name','Formula','P_Exact',
                  'P_Congruent','P_MinCon']
    metrics_cols = ['ClusterID','ClusterSize','Input','Name','Formula','Comps',
                    'N_Exact','N_Congruent','N_MinCon','P_Exact','P_Exact_sub',
                    'P_Exact_all','P_Congruent','P_Congruent_sub','P_Congruent_all',
                    'P_MinCon','P_MinCon_sub','P_MinCon_all']              

    """No arguments"""
    with pytest.raises(TypeError):
        SeroClust()
                
    """Invalid input - serovar parameter is empty list"""
    with pytest.raises(InvalidInput):
        SeroClust('clust1', [])

    """Invalid input - serovar parameter is not a list"""
    with pytest.raises(InvalidInput):
        SeroClust('clust1', '1')
        
    """Only invalid serovar"""
    with pytest.raises(InvalidInput):
        SeroClust('clust',WKLMSerovar('Test'))
               
    """Valid serovar name"""
    wklm_objs = [WKLMSerovar('Norwich')]
    assert_frame_equal(SeroClust('clust1',wklm_objs).results.reset_index(drop=True),
        pd.DataFrame({'ClusterID': ['clust1'], 'ClusterSize': [1], 'Input': ['Norwich'], 
                      'Name': ['Norwich'], 'Formula': ['I 6,7:e,h:1,6'],'P_Exact': [1.0], 
                      'P_Congruent': [1.0], 'P_MinCon': [1.0]},columns = results_cols).reset_index(drop=True))

    """Default - results"""
    wklm_objs = [WKLMSerovar(i) for i in ['Kivu','Kivu','Javiana']]
    assert_frame_equal(SeroClust('clust1',wklm_objs).results.reset_index(drop=True),
        pd.DataFrame({'ClusterID': ['clust1'], 'ClusterSize': [3], 'Input': ['Kivu'], 
                      'Name': ['Kivu'], 'Formula': ['I 6,7:d:1,6'],'P_Exact': [0.6667], 
                      'P_Congruent': [0.6667], 'P_MinCon': [0.6667]},columns = results_cols).reset_index(drop=True))

    """Default - metrics"""
    wklm_objs = [WKLMSerovar(i) for i in ['Kivu','Kivu','Javiana']]
    assert_frame_equal(SeroClust('clust1',wklm_objs).metrics.reset_index(drop=True),
        pd.DataFrame({'ClusterID': ['clust1','clust1'], 'ClusterSize': [3,3], 
                      'Input': ['Kivu','Javiana'], 'Name': ['Kivu','Javiana'], 
                      'Formula': ['I 6,7:d:1,6','I [1],9,12:l,z28:1,5:[R1…]'],
                      'Comps': [['exact', 'incongruent'],['incongruent', 'exact']],
                      'N_Exact':[2,1],'N_Congruent':[2,1],'N_MinCon':[2,1],
                      'P_Exact': [0.6667,0.3333],'P_Exact_sub': [0.6667,0.3333], 
                      'P_Exact_all': [0.6667,0.3333],'P_Congruent': [0.6667,0.3333],
                      'P_Congruent_sub': [0.6667,0.3333], 'P_Congruent_all': [0.6667,0.3333],
                      'P_MinCon': [0.6667,0.3333],'P_MinCon_sub': [0.6667,0.3333], 
                      'P_MinCon_all': [0.6667,0.3333]},columns = metrics_cols).reset_index(drop=True))
    
    """Default with formula"""
    wklm_objs = [WKLMSerovar(i) for i in ['Kivu','I 6,7:d:1,6','Javiana']]
    assert_frame_equal(SeroClust('clust1',wklm_objs).results.reset_index(drop=True),
        pd.DataFrame({'ClusterID': ['clust1'], 'ClusterSize': [3], 'Input': ['Kivu'], 
                      'Name': ['Kivu'], 'Formula': ['I 6,7:d:1,6'],
                      'P_Exact': [0.6667], 'P_Congruent': [0.6667], 
                      'P_MinCon': [0.6667]},columns = results_cols).reset_index(drop=True))
                 
    """Handle acceptable missing serovar notations""" 
    wklm_objs = [WKLMSerovar(i) for i in ['Austin','NULL','NULL','NA','']]
    assert_frame_equal(SeroClust('clust1',wklm_objs).results.reset_index(drop=True),
        pd.DataFrame({'ClusterID': ['clust1'], 'ClusterSize': [5], 'Input': ['Austin'], 
                      'Name': ['Austin'], 'Formula': ['I 6,7:a:1,7'],
                      'P_Exact': ['1.0(0.2)'], 'P_Congruent': ['1.0(0.2)'], 
                      'P_MinCon': ['1.0(0.2)']},columns = results_cols).reset_index(drop=True))

    """It's a tie""" 
    wklm_objs = [WKLMSerovar(i) for i in ['Austin','Kivu']]
    assert_frame_equal(SeroClust('clust1',wklm_objs).results.reset_index(drop=True),
        pd.DataFrame({'ClusterID': ['clust1','clust1'], 'ClusterSize': [2,2], 
                      'Input': ['Austin','Kivu'], 'Name': ['Austin','Kivu'], 
                      'Formula': ['I 6,7:a:1,7','I 6,7:d:1,6'],
                      'P_Exact': [0.5,0.5], 'P_Congruent': [0.5,0.5], 
                      'P_MinCon': [0.5,0.5]},columns = results_cols).reset_index(drop=True))
                 
    """It's a tie broken with a congruent serovar"""  
    wklm_objs = [WKLMSerovar(i) for i in ['Coeln','Coeln','I 1,4,12:y:1,2','Kivu','Kivu']]
    assert_frame_equal(SeroClust('clust1',wklm_objs).results.reset_index(drop=True),
        pd.DataFrame({'ClusterID': ['clust1'], 'ClusterSize': [5], 'Input': ['Coeln'], 
                      'Name': ['Coeln'], 'Formula': ['I [1],4,[5],12:y:1,2'],
                      'P_Exact': [0.4], 'P_Congruent': [0.6], 
                      'P_MinCon': [0.6]},columns = results_cols).reset_index(drop=True))

    """It's a tie broken by minimally congruent serovar"""
    sort_by = ['m','c','e'] 
    wklm_objs = [WKLMSerovar(i) for i in ['Coeln','Coeln','I 1,4,12:y:1','Kivu','Kivu']]
    assert_frame_equal(SeroClust('clust1',wklm_objs,sort_by=sort_by).results.reset_index(drop=True),
        pd.DataFrame({'ClusterID': ['clust1'], 'ClusterSize': [5], 'Input': ['Coeln'], 
                      'Name': ['Coeln'], 'Formula': ['I [1],4,[5],12:y:1,2'],
                      'P_Exact': [0.4], 'P_Congruent': [0.4], 
                      'P_MinCon': [0.6]},columns = results_cols).reset_index(drop=True))
    
    """Testing congruency calculations"""
    sort_by = ['m','c','e'] 
    wklm_objs = [WKLMSerovar(i) for i in ['Montevideo','Montevideo','Montevideo',
              'I 6,7,[14],[54]:g,m,[p],s:[1,2,7]','I 6,7,[14],[54]:g,m,[p],s:[1,2,7]',
              'I 6,7:g,m,[p],s:-','I 6,7:g,m,[p],s:-','I 6:g,m,[p],s:-','Javiana','NULL']]                
    assert_frame_equal(SeroClust('cluster',wklm_objs,sort_by=sort_by).results.reset_index(drop=True),
        pd.DataFrame({'ClusterID': ['cluster'], 'ClusterSize': [10], 'Input': ['Montevideo'], 
                      'Name': ['Montevideo'], 'Formula': ['I 6,7,[14],[54]:g,m,[p],s:[1,2,7]'],
                      'P_Exact': ['0.5556(0.5)'], 'P_Congruent': ['0.7778(0.7)'], 
                      'P_MinCon': ['0.8889(0.8)']},columns = results_cols).reset_index(drop=True))


def test_all_antigens_missing():

    assert st.all_antigens_missing(WKLMSerovar('III'))      == True
    assert st.all_antigens_missing(WKLMSerovar('i'))        == True
    assert st.all_antigens_missing(WKLMSerovar('I :'))      == True
    assert st.all_antigens_missing(WKLMSerovar('II ::'))    == True
    assert st.all_antigens_missing(WKLMSerovar('v –:–:–'))  == True
    assert st.all_antigens_missing(WKLMSerovar('test'))     == True
    assert st.all_antigens_missing(WKLMSerovar('I 9:'))     == False
    assert st.all_antigens_missing(WKLMSerovar('iv –:2,3')) == False
    
def test_cluster(tmpdir,capsys):                        

    """Input file - two clusters"""    
    file1 = tmpdir.join('sero_cluster.tsv')
    with open(str(file1), 'w') as f1:
        f1.write('{}\t{}\n{}\t{}\n'.format('clust1','Javiana',
                                           'clust2','Kumasi'))

    """Print results"""    
    st.cluster(input_file=file1)
    in_cap = capsys.readouterr()
    in_expected = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                  'ClusterID','ClusterSize','Input',  'Name',   'Formula',                   'P_Exact','P_Congruent','P_MinCon',
                  'clust1',    1,           'Javiana','Javiana','I [1],9,12:l,z28:1,5:[R1…]', 1.0,      1.0,          1.0, 
                  'clust2',    1,           'Kumasi', 'Kumasi', 'I 30:z10:e,n,z15',           1.0,      1.0,          1.0)
    assert in_cap.out == in_expected

    """Print results (v=2)"""    
    st.cluster(input_file=file1,v=2)
    in_cap2 = capsys.readouterr()
    in_expected2 = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                   'ClusterID','ClusterSize','Input',  'Name',   'Formula',                   'P_Exact','P_Congruent','P_MinCon',
                   'clust1',    1,            'Javiana','Javiana','I [1],9,12:l,z28:1,5:[R1…]',1.0,      1.0,          1.0, 
                   'clust2',    1,            'Kumasi', 'Kumasi', 'I 30:z10:e,n,z15',          1.0,      1.0,          1.0)
    assert in_cap2.out == in_expected2

    """Print serovars"""    
    st.cluster(input_file=file1,v=1)
    in_cap3 = capsys.readouterr()
    in_expected3 = '{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\n'.format(
                   'ClusterID','ClusterSize','Input',  'Name',   'Formula',                    
                   'clust1',    1,           'Javiana','Javiana','I [1],9,12:l,z28:1,5:[R1…]',
                   'clust2',    1,           'Kumasi', 'Kumasi', 'I 30:z10:e,n,z15')
    assert in_cap3.out == in_expected3
                                           
    """Print metrics"""    
    st.cluster(input_file=file1,v=3)
    in_cap4 = capsys.readouterr()
    in_expected4 = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                   'ClusterID','ClusterSize','Input',  'Name',   'Formula',                   'N_Exact','N_Congruent','N_MinCon','P_Exact','P_Exact_sub','P_Exact_all','P_Congruent','P_Congruent_sub','P_Congruent_all','P_MinCon','P_MinCon_sub','P_MinCon_all',
                   'clust1',    1,           'Javiana','Javiana','I [1],9,12:l,z28:1,5:[R1…]', 1,        1,            1,          1.0,      1.0,          1.0,          1.0,          1.0,              1.0,              1.0,       1.0,           1.0,
                   'clust2',    1,           'Kumasi', 'Kumasi', 'I 30:z10:e,n,z15',           1,        1,            1,          1.0,      1.0,          1.0,          1.0,          1.0,              1.0,              1.0,       1.0,           1.0)
    assert in_cap4.out == in_expected4
                           

def test_compare(tmpdir,capsys):                        

    """Input file"""    
    file1 = tmpdir.join('wklm_compare.tsv')
    with open(str(file1), 'w') as f1:
        f1.write('{}\t{}\n'.format('Kumasi','I 30:z10:e,n,z15'))
    
    st.compare(input_file=file1)
    in_cap = capsys.readouterr()
    in_expected = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                  'Serovar1','Name','Formula','Serovar2','Name','Formula','Result',
                  'Kumasi','Kumasi','I 30:z10:e,n,z15','I 30:z10:e,n,z15','Kumasi','I 30:z10:e,n,z15','exact')      
    assert in_cap.out == in_expected

    """Input file - with header"""
    file2 = tmpdir.join('wklm_compare_header.tsv')
    with open(str(file2), 'w') as f2:
        f2.write('{}\t{}\n'.format('Header1','Header2'))
        f2.write('{}\t{}\n'.format('Kumasi','I 30:z10:e,n,z15'))

    st.compare(input_file=file2, header=True)
    h_in_cap = capsys.readouterr()
    h_in_expected = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                  'Header1','Name','Formula','Header2','Name','Formula','Result',
                  'Kumasi','Kumasi','I 30:z10:e,n,z15','I 30:z10:e,n,z15','Kumasi','I 30:z10:e,n,z15','exact')      
    assert h_in_cap.out == h_in_expected
    
    """Command line args"""
    st.compare(subj='Kumasi',query='I 30:z10:e,n,z15')
    args_cap = capsys.readouterr()    
    args_expected = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                  'Serovar1','Name','Formula','Serovar2','Name','Formula','Result',
                  'Kumasi','Kumasi','I 30:z10:e,n,z15','I 30:z10:e,n,z15','Kumasi','I 30:z10:e,n,z15','exact')
    assert args_cap.out == args_expected
 
    """Invalid input""" 
    st.compare(subj='test',query='I 30:z10:e,n,z15')
    inv_cap = capsys.readouterr()    
    inv_expected = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                  'Serovar1','Name','Formula','Serovar2','Name','Formula','Result',
                  'test','NA','NA','I 30:z10:e,n,z15','Kumasi','I 30:z10:e,n,z15','invalid input')
    assert inv_cap.out == inv_expected
          

def test_fields_to_formula():

    """Default fields"""
    formula = ['I','3,{10}{15}{15,34}','z10','1,5','[z49]']
    assert st.fields_to_formula(formula) == 'I 3,{10}{15}{15,34}:z10:1,5:[z49]'

    """Nan subsp"""
    formula = [np.nan,'3,{10}{15}{15,34}','z10','1,5','[z49]'] 
    assert st.fields_to_formula(formula) == '3,{10}{15}{15,34}:z10:1,5:[z49]'

    """Empty subsp2"""
    formula = ['','3,{10}{15}{15,34}','z10','1,5','[z49]'] 
    assert st.fields_to_formula(formula) == '3,{10}{15}{15,34}:z10:1,5:[z49]'
    
    """Missing subsp"""
    formula = ['3,{10}{15}{15,34}','z10','1,5','[z49]']
    with pytest.raises(ValueError): 
        st.fields_to_formula(formula)

    """Nan other_H"""
    formula = ['I','3,{10}{15}{15,34}','z10','1,5',np.nan]
    assert st.fields_to_formula(formula) == 'I 3,{10}{15}{15,34}:z10:1,5'

    """Empty other_H"""
    formula = ['I','3,{10}{15}{15,34}','z10','1,5','']
    assert st.fields_to_formula(formula) == 'I 3,{10}{15}{15,34}:z10:1,5'

    """Missing other_H """
    formula = ['I','3,{10}{15}{15,34}','z10','1,5']
    with pytest.raises(ValueError):
        st.fields_to_formula(formula)

    """Nan antigens"""
    formula = ['I', np.nan, np.nan, np.nan, np.nan]
    assert st.fields_to_formula(formula) ==  'I –:–:–'

    """Missing antigens"""
    formula = ['I', '–', '–', '–', '']
    assert st.fields_to_formula(formula) ==  'I –:–:–'
        
    """Subsp only"""
    formula = ['I','','','','']
    assert st.fields_to_formula(formula) == 'I –:–:–'
       
    """Invalid input """
    with pytest.raises(ValueError):
        st.fields_to_formula('Enteritidis')


def test_find_matches(): 
 
    """Empty"""
    assert st.find_matches(WKLMSerovar('')) == []

    """Space"""
    assert st.find_matches(WKLMSerovar(' ')) == []

    """No subspecies, missing antigen"""
    assert st.find_matches(WKLMSerovar('–')) == []
    
#    """Subspecies only"""
#    assert len(st.find_matches(WKLMSerovar('I'))) == 1623
    
#     """All antigens missing"""
#    assert len(st.find_matches(WKLMSerovar('I :'))) == 1623
   
    """Query is a subset - missing antigen - all minimally congruent"""
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('I 1,9,12:b:–'))] \
        == ['I [1],9,12:b:1,2', 'I [1],9,12:b:1,5']

    """Query is a subset - missing antigen and optional factor - all minimally congruent"""
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('I [1],9,12:b:–'))] \
        == ['I 9,12:b:e,n,z15', 'I [1],9,12:b:1,2', 'I [1],9,12:b:1,5', 'I 9,12:b:1,7']

    """Query is a subset - missing antigen and optional factor removed - all minimally congruent"""
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('I 9,12:b:–'))] \
        == ['I 9,12:b:e,n,z15', 'I [1],9,12:b:1,2', 'I [1],9,12:b:1,5', 'I 9,12:b:1,7']

    """Testing results"""
    assert [obj.result for obj in st.find_matches(WKLMSerovar('I [1],4,12,27:r,[i]:e,n,z15'))] \
        == ['minimally congruent', 'exact']   

    """Query is a subset - with optional - exact and minimally congruent"""        
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('I 8,[20]:b:e,n,z15'))] \
        == ['I 8,[20]:b:e,n,z15', 'I 6,8:b:e,n,z15'] 
        
    """Query is a subset - all minimally congruent"""
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('I 6:d:1,5:[z58]'))] \
        == ['I 6,7,[14]:d:1,5', 'I 6,8:d:1,5:[z58]', 'I 6,14,24:d:1,5', 'I [1],6,14,[25]:d:1,5'] 

    """Query is a subset - congruent and minimally congruent"""
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('I 8:d:1,5:[z58]'))] \
        == ['I 8,[20]:d:1,5', 'I 6,8:d:1,5:[z58]'] 
 
    """Query is a superset - minimally congruent"""
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('I 6,8,10:d:1,5:[z58]'))] \
        == ['I 6,8:d:1,5:[z58]']
       
    """Missing subspecies"""
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('[1],4,[5],12,[27]:b:1,5'))] \
        == ['II 4,12:b:1,5', 'I [1],4,[5],12,[27]:b:1,5']   

    """Missing subspecies2"""
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('1,4,[5],12,[27]:b:1,5'))] \
	    == ['I [1],4,[5],12,[27]:b:1,5']   
        
    """Missing subspecies and antigen"""
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('9,12:b:–'))] \
        == ['I [1],9,12:b:1,2', 'I [1],9,12:b:1,5', 'I 9,12:b:1,7', 'II [1],9,12:b:e,n,x', 'I 9,12:b:e,n,z15', 'II [1],9,12:b:z6', 'II [1],9,12:b:z39', 'II 1,9,12,46,27:b:z39']

    """No matches"""        
    assert [obj.query.formula for obj in st.find_matches(WKLMSerovar('I 1,9,12:b:6'))] \
        == []

    """Variant"""        
    assert [obj.query.name for obj in st.find_matches(WKLMSerovar('I 6,7,[14]:l,v:z6'))] \
        == ['Gdansk', 'Gdansk var. 14+']
           
    """Serovars with identical formulas"""        
    assert [obj.query.name for obj in st.find_matches(WKLMSerovar('I 1,4,[5],12:b:1,2:[z5],[z33]'))] \
        == ['Paratyphi B var. L(+) tartrate (= d–tartrate)+', 'Paratyphi B']

    assert [obj.query.name for obj in st.find_matches(WKLMSerovar('Miami'))] \
        == ['Miami or Sendai']
     
    assert [obj.query.name for obj in st.find_matches(WKLMSerovar('Choleraesuis'))] \
        == ['Choleraesuis or Typhisuis', 'Paratyphi C']
    

def test_formula_to_fields():

    """Default formula"""
    formula = 'I 3,{10}{15}{15,34}:z10:1,5:[z49]'
    assert st.formula_to_fields(formula) == ['I','3,{10}{15}{15,34}','z10','1,5','[z49]']

    """No subsp"""
    formula = '3,{10}{15}{15,34}:z10:1,5:[z49]'
    assert st.formula_to_fields(formula) == [np.nan,'3,{10}{15}{15,34}','z10','1,5','[z49]']

    """Missing antigens"""
    formula = 'I ::'
    assert st.formula_to_fields(formula) == ['I','–','–','–','']
 
    """Missing antigens3"""
    formula = 'I –:–:–'
    assert st.formula_to_fields(formula) == ['I','–','–','–','']
       
    """Subsp only"""
    formula = 'I'
    assert st.formula_to_fields(formula) == ['I','–','–','–','']
       
    """Invalid formula """
    with pytest.raises(InvalidInput):
        st.formula_to_fields('Enteritidis')


def test_get_factor():

    """Not optional"""
    assert st.get_factor('2','1,2,3') == '2'

    """Optional"""
    assert st.get_factor('s','f,g,[s]') == '[s]'

    """Exclusive"""
    assert st.get_factor('g','f,{g},s') == '{g}'

    """Weakly agglutinable"""
    assert st.get_factor('f','(f),g,s') == '(f)'

    """Other"""
    assert st.get_factor('R1','[R1…],[z37],[z45],[z49]') == '[R1…]'

    
def test_input_to_wklm():

    """Serovar name"""
    assert st.input_to_wklm('Montevideo').formula == 'I 6,7,[14],[54]:g,m,[p],s:[1,2,7]'

    """Multiple serovars"""
    assert st.input_to_wklm('Miami or Sendai').formula == 'I [1],9,12:a:1,5'


def test_is_min_subset():

    assert st.is_min_subset('5','5') == True
    assert st.is_min_subset('5','4') == False
    
    """z15 is a subset of e,n,z15, e,n,z15 is not a subset of z15"""    
    assert st.is_min_subset('e,n,z15','z15') == False
    
    """Not a subset"""
    assert st.is_min_subset('e,n,z15','z39') == False
    
    """e,n,z15 is a proper subset of e,n,z15,z39"""
    assert st.is_min_subset('e,n,z15','e,n,z15,z39') == True
    
    """Not a subset"""
    assert st.is_min_subset('e,n,z15','e,n,z39') == False
    
    """Both are proper subsets"""
    assert st.is_min_subset('6,7,8,[14],[54]','6,7,8') == True
    
    """6,7,8 is a proper subset of 6,7,8,9"""
    assert st.is_min_subset('6,7,8,[14],[54]','6,7,8,9') == True
    
    """Not a subset"""    
    assert st.is_min_subset('6,7,8,10,[14],[54]','6,7,8,9') == False
    
    """Both are proper subsets"""
    assert st.is_min_subset('g,m,[p],s','g,m,p,s') == True
    
    """m,p,s is a proper subset of g,m,[p],s"""
    assert st.is_min_subset('m,p,s','g,m,[p],s') == True
 
    """g,m,[p],s is not a subset of m,p,s"""
    assert st.is_min_subset('g,m,[p],s','m,p,s') == False
   
    """m is a proper subset of g,[m],[s],[t]"""
    assert st.is_min_subset('m','g,[m],[s],[t]') == True
    
    """the empty set is a proper subset of the empty set"""
    assert st.is_min_subset('[1,2,7]','–') == True

    """the empty set is a proper subset of the empty set"""
    assert st.is_min_subset('–','[1,2,7]') == True
    
    """the empty set is a proper subset of 1,2,7"""
    assert st.is_min_subset('–','1,2,7') == True

    """the empty set is a proper subset of the empty set"""    
    assert st.is_min_subset('–','–') == True

    """the empty set is a proper subset of the empty set"""    
    assert st.is_min_subset('[1,2,7]','[5]') == True
    
    
def test_is_name():

    assert st.is_name('Enteritidis') == True
    assert st.is_name('Paratyphi B') == True
    assert st.is_name('Paratyphi B var. L(+) tartrate (= d–tartrate)+') == True
    assert st.is_name('Senftenberg or Dessau') == True
    assert st.is_name('Senftenberg/Dessau') == True
    assert st.is_name('IV Houten') == True
    assert st.is_name('I 1,2,12:a:[1,5]') == False
    assert st.is_name('I') == False
    assert st.is_name('I :') == False
    assert st.is_name('Paratyphi A (I 1,2,12:a:[1,5])') == False
    assert st.is_name('') == True   # no checks here
    assert st.is_name(' ') == True  # no checks here


def test_is_opt_factor():

    assert st.is_opt_factor('g','g,[m],(s),{t}') == False
    assert st.is_opt_factor('m','g,[m],(s),{t}') == True
    assert st.is_opt_factor('s','g,[m],(s),{t}') == True
    assert st.is_opt_factor('t','g,[m],(s),{t}') == True
    assert st.is_opt_factor('x','g,[m],(s),{t}') == False


def test_is_opt_subset():

    """Any differences must be due to optional factors"""    
    assert st.is_opt_subset('5','5') == True
    assert st.is_opt_subset('5','4') == False  
    assert st.is_opt_subset('e,n,z15','z15') == False
    assert st.is_opt_subset('e,n,z15','z39') == False
    assert st.is_opt_subset('e,n,z15','e,n,z15,z39') == False
    assert st.is_opt_subset('e,n,z15','e,n,z39') == False
    assert st.is_opt_subset('6,7,8,[14],[54]','6,7,8') == True
    assert st.is_opt_subset('6,7,8,[14],[54]','6,7,8,9') == False
    assert st.is_opt_subset('6,7,8,10,[14],[54]','6,7,8,9') == False
    assert st.is_opt_subset('g,m,[p],s','g,m,p,s') == True
    assert st.is_opt_subset('g,m,[p],s','g,m,s') == True
    assert st.is_opt_subset('g,m,[p],s','m,p,s') == False
    assert st.is_opt_subset('g,[m],[s],[t]','m') == False
    assert st.is_opt_subset('[1,2,7]','–') == True
    assert st.is_opt_subset('1,2,7','–') == False
    assert st.is_opt_subset('–','–') == True
    assert st.is_opt_subset('[1,2,7]','[5]') == True


def test_matching_indices():

    assert st.matching_indices('5',['1,2','2,5,6','10,15','2,5,7,10']) == [1,3]
    assert st.matching_indices('–',['1,2','2,5,6,z15','–','–']) == [2,3]
    assert st.matching_indices('i',['I','I','V','III']) == [0,1]

  
def test_max_factors():

    assert st.max_factors('[e,n,x]') == {'e', 'n', 'x'}
    assert st.max_factors('[e],(n),{x}') == {'e', 'n', 'x'}
    assert st.max_factors('1,2,5') == {'1', '2', '5'}
    assert st.max_factors('[1],2,5') == {'1', '2', '5'}
    assert st.max_factors('–') == {'–'}


def test_merge_wklm_objs(caplog):    
    
    sero_index=['Name','Std_Name','Formula','Std_Formula','Species','Subspecies',
            'O','P1','P2','other_H','Group','Old_Group','Input']
              
    """Serovar name"""
    assert st.input_to_wklm('Montevideo').formula == 'I 6,7,[14],[54]:g,m,[p],s:[1,2,7]'

    """Serovar formula"""
    assert st.input_to_wklm('I 6,7,[14],[54]:g,m,[p],s:[1,2,7]').name == 'Montevideo'
    
    """Serovars with identical formulas"""
    assert st.input_to_wklm('Miami or Sendai').meta.equals(pd.Series({ 
          'Name':'Miami or Sendai','Std_Name':'miami or sendai',
          'Formula':'I [1],9,12:a:1,5','Std_Formula':'i 1,9,12:a:1,5',
          'Species':'enterica','Subspecies':'I','O':'[1],9,12','P1':'a','P2':'1,5',
          'other_H':'','Group':'O:9','Old_Group':'D1','Input':'Miami or Sendai'},
          index=sero_index))
    
    """Invalid input - incongruent"""
    caplog.clear()
    assert st.input_to_wklm('I 6,7:y:1,6 or I 6,7:y:1,7').meta.equals(pd.Series({ 
          'Name':np.nan,'Std_Name':np.nan,
          'Formula':np.nan,'Std_Formula':np.nan,
          'Species':np.nan,'Subspecies':np.nan,'O':np.nan,'P1':np.nan,'P2':np.nan,
          'other_H':np.nan,'Group':np.nan,'Old_Group':np.nan,'Input':'I 6,7:y:1,6 or I 6,7:y:1,7'},
          index=sero_index))
    assert caplog.records[0].message == "The input serovars are incongruent and cannot be merged."

    """Merge minimally congruent"""
    assert st.input_to_wklm('I 6,7:y:1,6 or I 6,7:y:1').meta.equals(pd.Series({ 
          'Name':np.nan,'Std_Name':np.nan,
          'Formula':'I 6,7:y:1','Std_Formula':'i 6,7:y:1',
          'Species':'enterica','Subspecies':'I','O':'6,7','P1':'y','P2':'1',
          'other_H':'','Group':np.nan,'Old_Group':np.nan,'Input':'I 6,7:y:1,6 or I 6,7:y:1'},
          index=sero_index))
  
    """Merge minimally congruent - missing antigen"""
    assert st.input_to_wklm('I 6,7:l,w:1,5 or I 6,7:–:1,5').meta.equals(pd.Series({ 
          'Name':np.nan,'Std_Name':np.nan,
          'Formula':'I 6,7:–:1,5','Std_Formula':'i 6,7:–:1,5',
          'Species':'enterica','Subspecies':'I','O':'6,7','P1':'–','P2':'1,5',
          'other_H':'','Group':np.nan,'Old_Group':np.nan,'Input':'I 6,7:l,w:1,5 or I 6,7:–:1,5'},
          index=sero_index))
    
    """Different subsp, identical formulas"""
    caplog.clear()
    input = 'I [1],4,12,27:k:1,6 or II [1],4,12,27:k:1,6'
    assert st.input_to_wklm(input).meta.equals(pd.Series({ 
          'Name':np.nan,'Std_Name':np.nan,
          'Formula':np.nan,'Std_Formula':np.nan,
          'Species':np.nan,'Subspecies':np.nan,'O':np.nan,'P1':np.nan,'P2':np.nan,
          'other_H':np.nan,'Group':np.nan,'Old_Group':np.nan,'Input':'I [1],4,12,27:k:1,6 or II [1],4,12,27:k:1,6'},
          index=sero_index))
    assert caplog.records[0].message == "The input serovars are incongruent and cannot be merged."


def test_min_factors():

    assert st.min_factors('[e,n,x]') == set()
    assert st.min_factors('[e],(n),{x}') == set()
    assert st.min_factors('1,2,5') == {'1', '2', '5'}
    assert st.min_factors('[1],2,5') == {'2', '5'}
    assert st.min_factors('–') == {'–'}

 
def test_prep():

    assert st.prep('I 1,(4),[5],{10}{15}{15,34}:b:1,2:[z5],[z33]') == 'i 1,4,5,10,15,15,34:b:1,2:z5,z33'
    assert st.prep('I [1],4,[5],12:e,h:1,5:[R1…]') == 'i 1,4,5,12:e,h:1,5:r1'
    assert st.prep('I [1],4,[5],12:e,h:1,5:[R1...]') == 'i 1,4,5,12:e,h:1,5:r1'

    
def test_query(tmpdir,capsys):                        

    file = tmpdir.join('wklm_query.tsv')   
    with open(str(file), 'w') as f:
        f.write('{}\n'.format('Kumasi'))
    
    """Input file - single query""" 
    st.query(input_file=file)
    in_cap = capsys.readouterr()
    in_expected ='{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n'.format(
                 'Input','Name','Formula','Match',
                 'Kumasi','Kumasi','I 30:z10:e,n,z15','exact')       
    assert in_cap.out == in_expected

    file2 = tmpdir.join('wklm_query2.tsv')   
    with open(str(file2), 'w') as f2:
        f2.write('{}\n{}\n'.format('I 6,7,[14],[54]:g,m,[p],s:[1,2,7]','I [1],9,12:l,z28:1,5:[R1…]'))

    """Input file - multiple queries""" 
    st.query(input_file=file2)
    in_cap = capsys.readouterr()
    in_expected ='{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n'.format(
                 'Input','Name','Formula','Match',
                 'I 6,7,[14],[54]:g,m,[p],s:[1,2,7]','Montevideo','I 6,7,[14],[54]:g,m,[p],s:[1,2,7]','exact',
                 'I [1],9,12:l,z28:1,5:[R1…]','Javiana','I [1],9,12:l,z28:1,5:[R1…]','exact')
    assert in_cap.out == in_expected
    
    """Command line arg"""
    st.query(serovar='Kumasi')
    arg_cap = capsys.readouterr()    
    arg_expected ='{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n'.format(
                 'Input','Name','Formula','Match',
                 'Kumasi','Kumasi','I 30:z10:e,n,z15','exact',)
    assert arg_cap.out == arg_expected

    """Exact match - identical formulas, name"""
    st.query(serovar='Choleraesuis',exact=True)
    arg_cap = capsys.readouterr()    
    arg_expected ='{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n'.format(
                 'Input','Name','Formula','Match',
                 'Choleraesuis','Choleraesuis','I 6,7:c:1,5','exact',)
    assert arg_cap.out == arg_expected

    """Exact match - identical formulas, formula"""
    st.query(serovar='I 6,7:c:1,5',exact=True)
    arg_cap = capsys.readouterr()    
    arg_expected ='{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n'.format(
                 'Input','Name','Formula','Match',
                 'I 6,7:c:1,5','Choleraesuis or Typhisuis','I 6,7:c:1,5','exact',)
    assert arg_cap.out == arg_expected

    """Exact and congruent"""
    st.query(serovar='I 6,7:c:1,5')
    arg_cap = capsys.readouterr()    
    arg_expected ='{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n'.format(
                 'Input','Name','Formula','Match',
                 'I 6,7:c:1,5','Choleraesuis or Typhisuis','I 6,7:c:1,5','exact',
                 'I 6,7:c:1,5','Paratyphi C','I 6,7,[Vi]:c:1,5','congruent')
    assert arg_cap.out == arg_expected

    """No exact matches"""
    st.query(serovar='I 1,9,12:b:–',exact=True)
    arg_cap = capsys.readouterr()    
    arg_expected ='{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n'.format(
                 'Input','Name','Formula','Match',
                 'I 1,9,12:b:–','NA','NA','none')
    assert arg_cap.out == arg_expected

    """Only minimally congruent matches"""
    st.query(serovar='I 1,9,12:b:–')
    arg_cap = capsys.readouterr()    
    arg_expected ='{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n'.format(
                 'Input','Name','Formula','Match',
                 'I 1,9,12:b:–','Onarimon','I [1],9,12:b:1,2','minimally congruent',
                 'I 1,9,12:b:–','Frintrop','I [1],9,12:b:1,5','minimally congruent')
    assert arg_cap.out == arg_expected

    """Invalid input""" 
    st.query(serovar='test')
    inv_cap = capsys.readouterr()    
    inv_expected ='{}\t{}\t{}\t{}\n{}\t{}\t{}\t{}\n'.format(
                 'Input','Name','Formula','Match',
                 'test','NA','NA','none')
    assert inv_cap.out == inv_expected


def test_split_input():

    """Multiple serovar predictions as input"""
    assert st.split_input('Paratyphi C or Choleraesuis or Typhisuis*') \
        == ['Paratyphi C', 'Choleraesuis', 'Typhisuis*']
    assert st.split_input('Tananarive/Brunei') == ['Tananarive', 'Brunei']


def test_standardize_formula():

    """Missing antigens"""
    assert st.standardize_formula('I 4,5,12:nonmotile') == 'I 4,5,12:–:–'
    assert st.standardize_formula('I 4,5,12:Non-Motile') == 'I 4,5,12:–:–'
    assert st.standardize_formula('I Rough:i:1,2') == 'I –:i:1,2'
    assert st.standardize_formula('I Rough:non-motile') == 'I –:–:–'
    assert st.standardize_formula('I mucoid:i:1,2') == 'I –:i:1,2'
    assert st.standardize_formula('I 4,5,12::') == 'I 4,5,12:–:–'
    assert st.standardize_formula('I 4,5,12:i:undetermined') == 'I 4,5,12:i:–'
    assert st.standardize_formula('I Undetermined:i:1,2') == 'I –:i:1,2'
    assert st.standardize_formula('I :i:1,2') == 'I –:i:1,2'
    assert st.standardize_formula(':i:1,2') == '–:i:1,2'
    
    """Separate exclusive factors"""
    assert st.standardize_formula('I 1,4,5,{10}{15}{15,34}:b:1,2') == 'I 1,4,5,{10},{15},{15,34}:b:1,2'
    
    """Invalid input"""
    with pytest.raises(InvalidInput):
        st.standardize_formula('Rough:non-motile')
    with pytest.raises(InvalidInput):
        st.standardize_formula('–:–:–')
    with pytest.raises(InvalidInput):
        st.standardize_formula('–:–:–:–')
    

def test_standardize_input():

    assert st.standardize_input('Poona*') == 'Poona'

    """Remove content in parentheses, except for weakly agglutinable factors."""
    assert st.standardize_input('Cerro var. 14+ (Siegburg)') == 'Cerro var. 14+'
    assert st.standardize_input('I 1,3,10,19:f,g,t:1,(2),7') == 'I 1,3,10,19:f,g,t:1,(2),7'
    assert st.standardize_input('IIIb 16:(k):e,n,x,z15') == 'IIIb 16:(k):e,n,x,z15'
    
    """Transform subsp."""
    assert st.standardize_input('Enterica 16:i:z6') == 'I 16:i:z6'
    assert st.standardize_input('salamae 16:d:1,5') == 'II 16:d:1,5'
    assert st.standardize_input('Arizonae 13,23:z4,z24:–') == 'IIIa 13,23:z4,z24:–'
    assert st.standardize_input('Diarizonae (6),14:r:z') == 'IIIb (6),14:r:z'
    assert st.standardize_input('houtenae 6,14:z4,z23:–') == 'IV 6,14:z4,z23:–'
    assert st.standardize_input('bongori 1,40:z35:–') == 'V 1,40:z35:–'
    assert st.standardize_input('indica 41:b:1,7') == 'VI 41:b:1,7'
    
    """Remove leading/trailing whitespace"""
    assert st.standardize_input(' Enterica 16:i:z6 ') == 'I 16:i:z6'
    
 
def test_standardize_unicode():
 
    assert st.standardize_unicode('I [1],4,[5],12:e,h:1,5:[R1...]') == 'I [1],4,[5],12:e,h:1,5:[R1…]'
    assert st.standardize_unicode('I 4,5,12:i:-') == 'I 4,5,12:i:–'
