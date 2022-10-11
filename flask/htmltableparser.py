#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 19:53:19 2021

@author: naim
from https://srome.github.io/Parsing-HTML-Tables-in-Python-with-BeautifulSoup-and-pandas/
"""

import requests
import pandas as pd
from bs4 import BeautifulSoup

class HTMLTableParser:
    
     def parse_url(self, url):
         response = requests.get(url)
         soup = BeautifulSoup(response.text, 'lxml')
         return [(table['id'],self.parse_html_table(table))\
                 for table in soup.find_all('table')]  
 
     def parse_html_table(self, table):
         n_columns = 0
         n_rows=0
         column_names = []
 
         # Find number of rows and columns
         # we also find the column titles if we can
         for row in table.find_all('tr'):
             
             # Determine the number of rows in the table
             td_tags = row.find_all('td')
             if len(td_tags) > 0:
                 n_rows+=1
                 if n_columns == 0:
                     # Set the number of columns for our table
                     n_columns = len(td_tags)
                     
             # Handle column names if we find them
             th_tags = row.find_all('th') 
             if len(th_tags) > 0 and len(column_names) == 0:
                 for th in th_tags:
                     column_names.append(th.get_text())
 
         # Safeguard on Column Titles
         if len(column_names) > 0 and len(column_names) != n_columns:
             raise Exception("Column titles do not match the number of columns")
 
         columns = column_names if len(column_names) > 0 else range(0,n_columns)
         df = pd.DataFrame(columns = columns, index= range(0,n_rows))
         row_marker = 0
         for row in table.find_all('tr'):
             column_marker = 0
             columns = row.find_all('td')
             for column in columns:
                 df.iat[row_marker,column_marker] = column.get_text()
                 column_marker += 1
             if len(columns) > 0:
                 row_marker += 1
                 
         # Convert to float if possible
         for col in df:
             try:
                 df[col] = df[col].astype(float)
             except ValueError:
                 pass
         
         return df