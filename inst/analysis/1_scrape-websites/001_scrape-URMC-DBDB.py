#!/usr/bin/env python3

# Imports.
import os
from os.path import dirname
import requests
from bs4 import BeautifulSoup
import pandas as pd

# Directories.
here = os.getcwd()
root = dirname(here)
downdir = os.path.join(root,"downloads")

# Access the webpage.
URL = "https://www.dbdb.urmc.rochester.edu/associations/list"
session = requests.session()
page = session.get(URL,verify=False) #It is not recommended to set verify=False.
# But, this overcomes error when trying to reach out the website.

# Parse the webpage.
soup = BeautifulSoup(page.content,"lxml")

# Get data from table.
data = list()
table = soup.find_all('table')[1] # we are interested in second data table.
table_rows = table.find_all('tr') # tr = table rows
table_headers = [th.text for th in table_rows.pop(0).find_all('th')]
for row in table_rows:
    row_data = row.find_all('td') # td = table data
    row_values = [i.text for i in row_data]
    data.append(row_values)
# Ends loop to collect table data.

# Create a pandas df. 
df = pd.DataFrame(data)
df.columns = table_headers

# Save the data.
myfile = os.path.join(downdir,"rochester-dbdb-associations.csv")
df.to_csv(myfile)
