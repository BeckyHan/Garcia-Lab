import os

def mkdirs(path):
    """
    Generate folder (including the path) based on the input path

    Parameters
    ----------
    path : str
        The input path to be created
    """

    isExists = os.path.exists(path)

    if not isExists:
        os.makedirs(path)