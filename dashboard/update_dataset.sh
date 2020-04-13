#!/bin/sh

rm timeseries.json
wget https://pomber.github.io/covid19/timeseries.json

echo "let world_data = " > world_data.js
cat timeseries.json >> world_data.js

