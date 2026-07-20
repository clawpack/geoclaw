"""
Tests for netcdf_utils.compression_encoding.

Covers the mapping from the user-facing ``compression`` argument to an xarray
encoding fragment, including the ``chunksizes`` handling.
"""
import pytest

from clawpack.geoclaw.netcdf_utils import compression_encoding

pytestmark = [pytest.mark.python, pytest.mark.netcdf]


@pytest.mark.parametrize("value", [None, False, 0, {}])
def test_falsy_gives_no_compression(value):
    """None/False/0/empty-dict all mean 'no compression' -> empty fragment."""
    assert compression_encoding(value) == {}


def test_true_is_level1_shuffle():
    """True selects the recommended default: zlib level 1 + shuffle."""
    assert compression_encoding(True) == {
        "zlib": True, "complevel": 1, "shuffle": True,
    }


def test_int_selects_complevel():
    """An int selects that zlib complevel (with shuffle)."""
    assert compression_encoding(5) == {
        "zlib": True, "complevel": 5, "shuffle": True,
    }


def test_dict_passthrough():
    """A dict is used verbatim."""
    enc = {"zlib": True, "complevel": 4, "shuffle": False}
    assert compression_encoding(enc) == enc


def test_chunksizes_applied_only_when_compressing():
    """chunksizes is added when compression is enabled, and omitted otherwise."""
    assert "chunksizes" not in compression_encoding(None, chunksizes=(1, 2, 3))
    enc = compression_encoding(True, chunksizes=(1, 2, 3))
    assert enc["chunksizes"] == (1, 2, 3)


def test_dict_chunksizes_not_overridden():
    """A dict that already sets chunksizes keeps its own value."""
    enc = compression_encoding({"zlib": True, "chunksizes": (7, 8)},
                               chunksizes=(1, 2))
    assert enc["chunksizes"] == (7, 8)


def test_invalid_type_raises():
    """A nonsensical compression argument raises a clear error."""
    with pytest.raises(ValueError, match="compression must be"):
        compression_encoding("gzip")
