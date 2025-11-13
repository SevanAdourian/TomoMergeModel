# Test file to run unit tests on creation of Dataset
# and MergeMethods object and testing functionality

from Dataset import Dataset
from MergeMethods import MergeMethods
import pytest

# first series of tests: checking for invalid object creation

# bad filePath for Global Dataset
def test_badFilePath_1():
    with pytest.raises(ValueError):
        temp = Dataset(1, Dataset.GLOBAL)

def test_badFilePath_2():
    with pytest.raises(ValueError):
        temp = Dataset(None, Dataset.GLOBAL)

def test_badFilePath_3():
    with pytest.raises(ValueError):
        temp = Dataset("hello", Dataset.GLOBAL)

# bad model type for any Dataset
def test_badModelType_1():
    with pytest.raises(ValueError):
        temp = Dataset("Data\\netcdf\\semucb-2014-ucb-vs.nc", None)

def test_badModelType_2():
    with pytest.raises(ValueError):
        temp = Dataset("Data\\netcdf\\semucb-2014-ucb-vs.nc", 3, "Data\\netcdf\\semucb-2014-ucb-vs.nc")

# bad spline file path for regional datasets
def test_badSplinePath_1():
    with pytest.raises(ValueError):
        temp = Dataset("Data\\netcdf\\semucb-2014-ucb-vs.nc", Dataset.REGIONAL, 1)

def test_badSplinePath_2():
    with pytest.raises(ValueError):
        temp = Dataset("Data\\netcdf\\semucb-2014-ucb-vs.nc", Dataset.REGIONAL)

def test_badSplinePath_3():
    with pytest.raises(ValueError):
        temp = Dataset("Data\\netcdf\\semucb-2014-ucb-vs.nc", Dataset.REGIONAL, "hello")

