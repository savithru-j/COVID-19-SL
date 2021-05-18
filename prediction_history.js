// COVID-19 Simulation Tool in JavaScript
// Copyright 2020 Savithru Jayasinghe and Dushan Wadduwage
// Licensed under the MIT License (LICENSE.txt)

'use strict';

var date_format = 'YYYY-MM-DD';

// [+0 days, +1 day, +2 days, +3 days, +5 days, +7 days, +14 days]
var prediction_data = [
  {t: "2021-05-16", y: [2431, 2445, 2457, 2471, 2497, 2527, 2624]},
  {t: "2021-05-04", y: [1949, 2055, 2158, 2276, 2518, 2790, 3987]},
  {t: "2021-05-01", y: [1586, 1745, 1921, 2116, 2564, 3105, 6071]},
];

function populatePredictionHistoryTable()
{
  let table = document.getElementById("table_pred_hist");

  let actual_data = world_data["Sri Lanka"];
  let t_offsets = [0, 1, 2, 3, 5, 7, 14];

  for (let i = 0; i < prediction_data.length; ++i)
  {
    const data = prediction_data[i];
    let pred_date = moment(data.t, date_format);
    let ind = -1;
    for (let j = 0; j < actual_data.length; ++j)
      if (moment(actual_data[j].date, date_format).isSame(pred_date))
      {
        ind = j;
        break;
      }
    console.log("ind: " + ind);

    let actual_new_cases = [];
    for (let j = 0; j < t_offsets.length; ++j)
    {
      let tmp = actual_data[ind + t_offsets[j]];
      if (tmp != null)
        actual_new_cases.push(tmp.confirmed - actual_data[ind + t_offsets[j] - 1].confirmed); //diff with previous day
    }

    {
      let row = table.insertRow();
      let cell0 = row.insertCell(0);
      cell0.innerHTML = data.t;
      cell0.rowSpan = 3;

      let cell1 = row.insertCell(1);
      cell1.innerHTML = "Predicted new cases";

      for (let j = 0; j < data.y.length; ++j)
      {
        let cellj = row.insertCell(2+j);
        cellj.innerHTML = data.y[j];
      }
    }
    {
      let row = table.insertRow();
      let cell0 = row.insertCell(0);
      cell0.innerHTML = "Actual new cases";

      for (let j = 0; j < data.y.length; ++j)
      {
        let cellj = row.insertCell(1+j);
        if (j < actual_new_cases.length)
          cellj.innerHTML = actual_new_cases[j];
        else
          cellj.innerHTML = '-';
      }
    }
    {
      let row = table.insertRow();
      let cell0 = row.insertCell(0);
      cell0.innerHTML = "% Difference";

      for (let j = 0; j < data.y.length; ++j)
      {
        let cellj = row.insertCell(1+j);
        if (j < actual_new_cases.length)
          cellj.innerHTML = ((data.y[j] - actual_new_cases[j]) / actual_new_cases[j] * 100.0).toFixed(2) + "%";
        else
          cellj.innerHTML = '-';
      }
    }



    // console.log(data);
  }

  // let row = table.insertRow();
  // for (let i = 0; i < 9; ++i)
  // {
  //   let cell = row.insertCell(i);
  //   cell.innerHTML = i;
  //   console.log(cell);
  // }

}
