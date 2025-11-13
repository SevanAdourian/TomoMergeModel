# Test file to run unit tests on creation of Dataset
# and MergeMethods object and testing functionality

from Dataset import Dataset
from MergeMethods import MergeMethods
import pytest

# first series of tests: checking for invalid object creation

# bad filePath for Dataset
def badFilePath_1():
    temp = Dataset(1, Dataset.GLOBAL)

def badFilePath_2():
    temp = Dataset(None, Dataset.GLOBAL)

def badFilePath_3():
    temp = Dataset("hello", Dataset.GLOBAL)

def test_badFilePath():
    with pytest.raises(ValueError):
        badFilePath_1()
        badFilePath_2()
        badFilePath_3()
