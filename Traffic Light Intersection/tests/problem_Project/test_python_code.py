# -*- coding: utf-8 -*-
# ASSIGNMENT: Project 1
# PROBLEM NUMBER: Project

# place as problem_x/test_name.py so that relative imports work

import pytest

from ..tst import _test_fileregex

POINTS = 28

def test_python_code(pattern=r'Submission/.*py(nb)*'):
    return _test_fileregex(pattern)


