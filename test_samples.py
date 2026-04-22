# Test file to run unit tests on creation of Dataset
# and MergeMethods object and testing functionality

from geostitch.Dataset import Dataset
from geostitch.MergeMethods import MergeMethods
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


def test_badFilePath_4():
    with pytest.raises(ValueError):
        temp = Dataset("", Dataset.GLOBAL)


# bad model type for any Dataset
def test_badModelType_1():
    with pytest.raises(ValueError):
        temp = Dataset("Data\\netcdf\\semucb-2014-ucb-vs.nc", None)


def test_badModelType_2():
    with pytest.raises(ValueError):
        temp = Dataset(
            "Data\\netcdf\\semucb-2014-ucb-vs.nc",
            3,
            "Data\\netcdf\\semucb-2014-ucb-vs.nc",
        )


# bad spline file path for regional datasets
def test_badSplinePath_1():
    with pytest.raises(ValueError):
        temp = Dataset("Data\\netcdf\\semucb-2014-ucb-vs.nc", Dataset.REGIONAL, 1, "m")


def test_badSplinePath_2():
    with pytest.raises(ValueError):
        temp = Dataset("Data\\netcdf\\semucb-2014-ucb-vs.nc", Dataset.REGIONAL)


def test_badSplinePath_3():
    with pytest.raises(ValueError):
        temp = Dataset(
            "Data\\netcdf\\semucb-2014-ucb-vs.nc", Dataset.REGIONAL, "hello", "m"
        )


# bad depth units
def test_badDepthUnits_1():
    with pytest.raises(ValueError):
        temp = Dataset(
            "./data/semucb-2014-ucb-vs.nc",
            Dataset.REGIONAL,
            "./data/CANVAS_15-60s_400km.nc",
            depthUnits="m",
        )


# test instances to check that parsing was done correctly
test_globalMod = Dataset(
    "./data/semucb-2014-ucb-vs.nc", Dataset.GLOBAL, depthUnits="m"
)  # could not parse value name
test_regionalMod = Dataset(
    "./data/CANVAS_15-60s_400km.nc",
    Dataset.REGIONAL,
    depthUnits="m",
)


# test_globalMod.plot_all_variables()

# bad model inputs for Merge


# No models provided
def test_badModel_1():
    with pytest.raises(ValueError):
        temp = MergeMethods(None, None, "conf.yaml")


# Only one model provided
def test_badModel_2():
    with pytest.raises(ValueError):
        temp = MergeMethods(None, test_globalMod, "conf.yaml")


def test_badModel_3():
    with pytest.raises(ValueError):
        temp = MergeMethods(test_globalMod, None, "conf.yaml")


# two global models provided
def test_badModel_4():
    with pytest.raises(ValueError):
        temp = MergeMethods(test_globalMod, test_globalMod, "conf.yaml")


# two regional models provided
def test_badModel_5():
    with pytest.raises(ValueError):
        temp = MergeMethods(test_regionalMod, test_regionalMod, "conf.yaml")


# bad conf file input for Merge


# No conf file path given
def test_badConfPath_1():
    with pytest.raises(ValueError):
        temp = MergeMethods(test_regionalMod, test_globalMod, None)


# non-existent file path
def test_badConfPath_2():
    with pytest.raises(ValueError):
        temp = MergeMethods(test_regionalMod, test_globalMod, "non-existent.yaml")


# valid non-yaml file path
def test_badConfPath_3():
    with pytest.raises(ValueError):
        temp = MergeMethods(test_regionalMod, test_globalMod, "README.md")


# valid yaml file path, but not the right variables
def test_badConfPath_4():
    with pytest.raises(ValueError):
        temp = MergeMethods(test_regionalMod, test_globalMod, "environment.yaml")
