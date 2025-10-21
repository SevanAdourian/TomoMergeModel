# imports here
import os, yaml

class Dataset():
    REGIONAL = 0
    GLOBAL = 1
    # constructors
    def __init__(self, filePath: str, modelType: int):

        # Input validation for file path
        if filePath is None:
            raise ValueError("Need to pass in file path")
        if type(filePath) is not str:
            raise ValueError("File path must be a string")
        if not os.path.exists(filePath):
            raise ValueError("Incorrect file path: file could not be found")
        
        # Input validation for model type
        if modelType not in [Dataset.REGIONAL, Dataset.GLOBAL]:
            raise ValueError("Model type must be Dataset.REGIONAL or Dataset.GLOBAL")
        
        # assign the file path and model type, and parse file for other attributes
        self.filePath = filePath
        self.modelType = modelType
        # self.parseModel(self.filePath, self.modelType) <--- To be implemented

    # Initialize a Dataset from the conf.yaml file
    def initFromConf(modelType: int):

        # Input validation for model type
        if modelType not in [Dataset.REGIONAL, Dataset.GLOBAL]:
            raise ValueError("Model type must be Dataset.REGIONAL or Dataset.GLOBAL")
        
        # load the conf.yaml file
        try:
            with open("conf.yaml", "r") as file:
                conf = yaml.safe_load(file)
        except OSError:
            raise ValueError("conf.yaml not found or not in working directory")
        # return to the user, a Dataset instance depending on the model type and file path in conf.yaml
        if modelType == Dataset.REGIONAL:
            return Dataset(conf["path_to_regional_model"], modelType)
        if modelType == Dataset.GLOBAL:
            return Dataset(conf["path_to_global_model"], modelType)


    # getters and setters here
    def getFilePath(self) -> str:
        return self.filePath

    def setFilePath(self, path: str):
        if path is None:
            raise ValueError("Need to pass in file path")
        if not isinstance(path, str):
            raise ValueError("File path must be a string")
        if not os.path.exists(path):
            raise ValueError("Incorrect file path: file could not be found")
        self.filePath = path

    def getModelType(self) -> int:
        return self.modelType

    def setModelType(self, mt: int):
        if mt not in [Dataset.REGIONAL, Dataset.GLOBAL]:
            raise ValueError("Model type must be Dataset.REGIONAL or Dataset.GLOBAL")
        self.modelType = mt


    # dataset specific functions here
    pass