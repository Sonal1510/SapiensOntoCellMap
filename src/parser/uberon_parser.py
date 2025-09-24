#!/usr/bin/python3

"""
Author          : Sonal Rashmi (expert review by Gemini)
Date            : 15/08/2025
Description     : Parses Uberon database and get the developmental layer, organs and tissues
"""

import os
import sys
import pandas as pd
from typing import Optional, List, Union

# Add the project root to sys.path to resolve absolute imports from the config directory.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from config.config import (
    UBERON_OBO_FILENAME,
    RAW_DATA_DIR
)

class UberonDBParser:
    """
    A class to parse and provide lookup functionalities for Uberon anatomical database
    """
    def __init__(self):
        """
        Initializes the UberonDBParser by loading and processing the file.
        It constructs file paths and prepares lookup structures for efficient querying.
        """
        self.uberon_filename = os.path.join(RAW_DATA_DIR, UBERON_OBO_FILENAME)
        self.uberon_