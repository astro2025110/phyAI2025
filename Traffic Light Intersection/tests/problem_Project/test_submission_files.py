# -*- coding: utf-8 -*-
# ASSIGNMENT: Project 1
# PROBLEM NUMBER: Project

# place as problem_x/test_name.py so that relative imports work

import pytest

from ..tst import _test_fileregex

POINTS = 70

def test_submission_files(pattern=r'Submission/project_1\.(pdf|ipynb)'):
    return _test_fileregex(pattern)


