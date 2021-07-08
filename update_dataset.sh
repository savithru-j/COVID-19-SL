#!/bin/sh

rm timeseries.json
wget https://pomber.github.io/covid19/timeseries.json

echo "let world_data = " > world_data.js
cat timeseries.json >> world_data.js

rm timeseries.json

./update_vaccine_data.py
