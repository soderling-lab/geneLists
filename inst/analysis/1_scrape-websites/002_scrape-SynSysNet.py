#!/usr/bin/env python3
''' Scrape Synaptic Protein-Tree from SynSysNet Database'''

# Imports.
import os
import re
import requests
import pandas as pd
from os.path import dirname
from bs4 import BeautifulSoup

# Directories.
here = os.getcwd()
root = dirname(dirname(here))
downdir = os.path.join(root,"downloads")

# Access the webpage.
URL = "http://bioinformatics.charite.de/synsys/index.php?site=syn_class"
session = requests.session()
page = session.get(URL)

# Parse the webpage.
soup = BeautifulSoup(page.content,"lxml")

# Get links to protein webpages.
BASE_URL = "http://bioinformatics.charite.de/synsys/"
regex = re.compile("./index.php\?site=search")
hrefs = soup.find_all('a',attrs={'href':regex})

# For every link, surf its page and collect all uniprot ids as comma separated string.
results = list()
for link in hrefs:
    url = BASE_URL + link.get('href').lstrip("./")
    webpage = session.get(url)
    soup = BeautifulSoup(webpage.content,"lxml")
    regex = re.compile("http://www.uniprot.org/.")
    proteins = soup.find_all('a',attrs={'href':regex})
    uniprot = [prot.get('href').split('uniprot/')[1] for prot in proteins]
    uniprot = ",".join(uniprot)
    protein_class = url.split("synapt=")[1]
    protein_subclass = url.split("Function=")[1].split("&")[0]
    row_key = "{}:{}".format(protein_class,protein_subclass)
    results.append({row_key:uniprot})
# Done.

# Merge results into single dictionary.
protein_tree = {}
for d in results:
    protein_tree.update(d)

# Create a pandas df. 
df = pd.DataFrame.from_dict(protein_tree,orient="index")

# Save the data.
myfile = os.path.join(downdir,"SynSysNet-Synaptic-Protein-Tree.csv")
df.to_csv(myfile)
