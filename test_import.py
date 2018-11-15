#!/usr/bin/env python3

"""
Created on Wed Sep 19 13:33:08 2018

@author: Xiaolong Cao
"""

import os
import sys
import argparse

# add current folder to system path
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from utils import test
from utils import countGapEachPositionFromMapFiles

if __name__ == '__main__':
    test.f()