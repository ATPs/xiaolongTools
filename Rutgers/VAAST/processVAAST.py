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
        return individuals
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
    

def getAbnormalIndividual(result=[], positive=[], negative=[], neutral=[], neutralAsPositive=True):
    '''
    positive is a list of positive individuals, or a str
    negative is a list of negative individuals, or a str
    neutral is a list of netral individuals, or a str
    if neutralAsPositive, neutral is considered positive, otherwise as negative
    result is a list of observed positives
    print the description about the result, which is positive is missing in result and which negative is observed in results
    for str of individuals, the str should be like "0,2,5-8", which means [0,2,5,6,7,8]
    '''
    result = strIndi2list(result)
    positive = strIndi2list(positive)
    negative = strIndi2list(negative)
    
    # neutralAsPositive
    if neutralAsPositive:
        positive = positive + neutral
    else:
        negative = negative + neutral
    
    # false negative, positive missing in result
    falseNegative = [e for e in positive if e not in result]
    
    # false positive, negative that in result
    falsePositive = [e for e in negative if e in result]
    
    to_print = []
    if len(falseNegative) !=0:
        to_print.append('false negative' + str(falseNegative))
    if len(falsePositive) !=0:
        to_print.append('false positive' + str(falsePositive))
    print(','.join(to_print))
    
    return falseNegative, falsePositive