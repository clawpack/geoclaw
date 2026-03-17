"""Configuration for PyTest"""

import pytest

def pytest_addoption(parser):
    parser.addoption('--save', action='store_true')

@pytest.fixture
def save(request):
    return request.config.getoption("--save", False)
