// COVID-19 Simulation Tool in JavaScript
// Copyright 2020 Savithru Jayasinghe and Dushan Wadduwage
// Licensed under the MIT License (LICENSE.txt)

'use strict';

var date_format = 'YYYY-MM-DD';


function populatePredictionHistoryTable()
{
  let table = document.getElementById("table_pred_hist");

  let row = table.insertRow();
  for (let i = 0; i < 9; ++i)
  {
    let cell = row.insertCell(i);
    cell.innerHTML = i;
    console.log(cell);
  }

}
