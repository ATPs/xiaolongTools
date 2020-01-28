# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 14:58:41 2019

@author: ATPs
"""

def strIndi2list(individuals):
    '''
    individuals is a list or a str.
    if str of individuals, the str should be like "0,2,5-8", which means [0,2,5,6,7,8]
    return a list
    '''
    if isinstance(individuals, list):
        return [int(e) for e in individuals]
    if isinstance(individuals, str):
        l = individuals.split(',')
        result = []
        for e in l:
            if '-' not in e:
                result.append(int(e))
            else:
                start, end = e.split('-')
                start = int(start)
                end = int(end) + 1
                result += list(range(start, end))
        return result
    print('wrong format', individuals)
    return None
    

def processVAASTGenoType(genotype, positive=None, negative=None):
    '''
    positive is a list of positive individuals, or a str like "3-5,7-9,11-12"
    negative is a list of negative individuals, or a str
    genotype is a genotype output of VAAST/pVAAST. it looks like "N|6|A:A|*:*;B|3-5,7-9,11-12|A:G|*:Q" or "0,2|-:C|-:G;3-6|-:C|-:G;7-8|-:C|-:G;10,12|-:C|-:G;13|-:C|-:G"
    '''
    if positive is None: positive = []
    if negative is None: negative = []
    positive = strIndi2list(positive)
    negative = strIndi2list(negative)
    GTs = genotype.split(';')
    detected_positive = []
    unknown = []
    for GT in GTs:
        if GT.count('|') == 3:
            individuals = GT.split('|')[1]
            individuals = strIndi2list(individuals)
            g = GT.split('|')[2]
        if GT.count('|') == 2:
            individuals = GT.split('|')[0]
            g = GT.split('|')[1]
        individuals = strIndi2list(individuals)
        
        if g == '^|^':
            unknown += individuals
        else:
            detected_positive += individuals
        
    detected_negative = [e for e in positive + negative if e not in detected_positive + unknown]
    true_positive = [e for e in detected_positive if e in positive]
    false_positive = [e for e in detected_positive if e not in positive]
    true_negative = [e for e in detected_negative if e in negative]
    false_negative = [e for e in detected_negative if e not in negative]
    
    return true_positive, true_negative, false_positive, false_negative, unknown