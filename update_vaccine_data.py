#! /usr/bin/python3

import wget 
import csv
import json
import os

csv_file_name = "vaccine_data.csv";

url = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vSckZ9FalsiAE1J8XbRKMCX_BFVYlvGPUlbddAki12uUFNoqc8lTucv3rSGn8ARsnIcM7c2Y_o4MQYB/pub?output=csv'
wget.download(url, csv_file_name)

csvfile = open(csv_file_name, 'r')
jsonfile = open('world_vaccine_data.js', 'w')

fieldnames = ("Date","Covishield 1","Covishield 2","Sinopharm 1","Sinopharm 2","Sputnik 1","Sputnik 2","Pfizer 1","Pfizer 2")
csvreader = csv.DictReader( csvfile, fieldnames)

next(csvreader) # This skips the first row of the CSV file.

jsonfile.write('let world_vaccine_data = {"Sri Lanka": [\n')

for row in csvreader:
    row_sum = {'t': row['Date'], 
               'dose1': (int(row['Covishield 1']) + int(row['Sinopharm 1']) + int(row['Sputnik 1']) + int(row['Pfizer 1'])),
               'dose2': (int(row['Covishield 2']) + int(row['Sinopharm 2']) + int(row['Sputnik 2']) + int(row['Pfizer 2']))
              }
    json.dump(row_sum, jsonfile)
    jsonfile.write(',\n')
    
jsonfile.write(']\n}')

if os.path.exists(csv_file_name):
  os.remove(csv_file_name)

