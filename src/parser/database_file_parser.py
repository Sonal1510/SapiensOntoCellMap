#!/usr/bin/python3

"""
Author         : Sonal Rashmi (expert review by Gemini)
Date           : 15/08/2025
Description    : Manages the parsing of any file type of database.
"""

import pandas as pd
import os
import sys

class DatabaseFileParser:
    """
    Parses a data file upon instantiation and stores the result
    in a 'dataframe' attribute.
    """
    def __init__(self, db_file, content_type):
        """
        Initializes the object by parsing the specified file.
        The resulting DataFrame is stored in self.dataframe.
        """
        print(f"--- Creating parser and parsing {os.path.basename(db_file)} ---")
        
        # Initialize the attribute to None
        self.dataframe = None 

        try:
            df = None
            if content_type == 'xlsx':
                df = pd.read_excel(db_file)
            elif content_type == 'csv':
                df = pd.read_csv(db_file)
            elif content_type == 'tsv':
                df = pd.read_csv(db_file, sep='\t')
            elif content_type == 'xlsx_multi_sheet':
                # Returns a Dict of DataFrames {sheet_name: dataframe}
                df = pd.read_excel(db_file, sheet_name=None)
            elif content_type == 'json':
                try:
                    df = pd.read_json(db_file)
                except ValueError:
                    df = pd.read_json(db_file, lines=True)
            
            if df is not None:
                # FIX: Handle printing logic based on whether df is a Dict or DataFrame
                if isinstance(df, dict):
                    print(f"✅ Successfully parsed multi-sheet file. Sheets found: {len(df.keys())}")
                else:
                    print(f"✅ Successfully parsed. Shape: {df.shape}")
                
                # Store the result in the instance attribute
                self.dataframe = df 
            else:
                print(f"Parser not implemented for content type '{content_type}'.")

        except Exception as e:
            print(f"❌ An error occurred while parsing: {e}")
            # self.dataframe remains None on failure