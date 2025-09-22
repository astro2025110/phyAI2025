# -*- coding: utf-8 -*-
# ASSIGNMENT: Project 2
# PROBLEM NUMBER: Project

# place as problem_x/test_name.py so that relative imports work

import pytest

from ..tst import _test_fileregex

POINTS = 2

def test_README(pattern=r'Submission/README.*'):
    return _test_fileregex(pattern)


