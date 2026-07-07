"""Configuration for PyTest"""

import pytest

def pytest_addoption(parser):
    parser.addoption('--save', action='store_true')

@pytest.fixture
def save(request):
    return request.config.getoption("--save", False)

@pytest.fixture(scope="session")
def download_cache(tmp_path_factory):
    r"""Session-scoped directory for files downloaded during testing.

    Tests that need a remote file (e.g. topography) should download it here
    rather than into a test-local ``tmp_path`` or the shared, persistent
    ``$CLAW/geoclaw/scratch`` directory used by non-test workflows.  Sharing
    this directory across tests in a session avoids re-downloading the same
    large reference files repeatedly, while still keeping it isolated to a
    single pytest session so it can never go stale or race with other runs.
    """
    return tmp_path_factory.mktemp("geoclaw_downloads")
