// COVID-19 Simulation Tool in JavaScript
// Copyright 2020 Savithru Jayasinghe and Dushan Wadduwage
// Licensed under the MIT License (LICENSE.txt)

'use strict';

var date_format = 'DD-MM-YYYY';

function getRawDataSriLanka()
{
  //Data array: [in_ward, recovered, deaths, foreign_input_to_quarantine]
  let data = [];
  data.push({t: moment('01-03-2020', date_format), y: [0, 1, 0, 0]});
  data.push({t: moment('02-03-2020', date_format), y: [0, 1, 0, 0]});
  data.push({t: moment('03-03-2020', date_format), y: [0, 1, 0, 0]});
  data.push({t: moment('04-03-2020', date_format), y: [0, 1, 0, 0]});
  data.push({t: moment('05-03-2020', date_format), y: [0, 1, 0, 0]});
  data.push({t: moment('06-03-2020', date_format), y: [0, 1, 0, 0]});
  data.push({t: moment('07-03-2020', date_format), y: [0, 1, 0, 0]});
  data.push({t: moment('08-03-2020', date_format), y: [0, 1, 0, 0]});
  data.push({t: moment('09-03-2020', date_format), y: [0, 1, 0, 0]});
  data.push({t: moment('10-03-2020', date_format), y: [0, 1, 0, 0]});
  data.push({t: moment('11-03-2020', date_format), y: [1, 1, 0, 0]});
  data.push({t: moment('12-03-2020', date_format), y: [2, 1, 0, 0]});
  data.push({t: moment('13-03-2020', date_format), y: [5, 1, 0, 2]});
  data.push({t: moment('14-03-2020', date_format), y: [10, 1, 0, 2]});
  data.push({t: moment('15-03-2020', date_format), y: [18, 1, 0, 7]});
  data.push({t: moment('16-03-2020', date_format), y: [28, 1, 0, 6]});
  data.push({t: moment('17-03-2020', date_format), y: [41, 1, 0, 4]});
  data.push({t: moment('18-03-2020', date_format), y: [52, 1, 0, 1]});
  data.push({t: moment('19-03-2020', date_format), y: [65, 1, 0, 2]});
  data.push({t: moment('20-03-2020', date_format), y: [71, 1, 0, 3]});
  data.push({t: moment('21-03-2020', date_format), y: [77, 1, 0, 6]});
  data.push({t: moment('22-03-2020', date_format), y: [86, 1, 0, 0]});
  data.push({t: moment('23-03-2020', date_format), y: [95, 2, 0, 0]});
  data.push({t: moment('24-03-2020', date_format), y: [99, 3, 0, 0]});
  data.push({t: moment('25-03-2020', date_format), y: [99, 3, 0, 0]});
  data.push({t: moment('26-03-2020', date_format), y: [99, 7, 0, 0]});
  data.push({t: moment('27-03-2020', date_format), y: [97, 9, 0, 0]});
  data.push({t: moment('28-03-2020', date_format), y: [103, 9, 1, 0]});
  data.push({t: moment('29-03-2020', date_format), y: [105, 11, 1, 0]});
  data.push({t: moment('30-03-2020', date_format), y: [106, 14, 2, 0]});
  data.push({t: moment('31-03-2020', date_format), y: [124, 17, 2, 0]});
  return data;
}

function getDataSriLanka()
{
  let data = [];
  for (let i = 0; i < data_raw_SL.length; ++i)
  {
    let num_confirmed_cases = 0;
    for (let j = 0; j < 3; ++j)
      num_confirmed_cases += data_raw_SL[i].y[j];
    data.push({t: data_raw_SL[i].t, y: num_confirmed_cases});
  }
  return data;
}

function getDataItaly()
{
  let data = [];
  data.push({t: moment('31-01-2020', date_format), y: 2});
  data.push({t: moment('06-02-2020', date_format), y: 3});
  data.push({t: moment('21-02-2020', date_format), y: 20});
  data.push({t: moment('22-02-2020', date_format), y: 79});
  data.push({t: moment('23-02-2020', date_format), y: 150});
  data.push({t: moment('24-02-2020', date_format), y: 227});
  data.push({t: moment('25-02-2020', date_format), y: 320});
  data.push({t: moment('26-02-2020', date_format), y: 445});
  data.push({t: moment('27-02-2020', date_format), y: 650});
  data.push({t: moment('28-02-2020', date_format), y: 888});
  data.push({t: moment('29-02-2020', date_format), y: 1128});
  data.push({t: moment('01-03-2020', date_format), y: 1694});
  data.push({t: moment('02-03-2020', date_format), y: 2036});
  data.push({t: moment('03-03-2020', date_format), y: 2502});
  data.push({t: moment('04-03-2020', date_format), y: 3089});
  data.push({t: moment('05-03-2020', date_format), y: 3858});
  data.push({t: moment('06-03-2020', date_format), y: 4636});
  data.push({t: moment('07-03-2020', date_format), y: 5883});
  data.push({t: moment('08-03-2020', date_format), y: 7375});
  data.push({t: moment('09-03-2020', date_format), y: 9172});
  data.push({t: moment('10-03-2020', date_format), y: 10149});
  data.push({t: moment('11-03-2020', date_format), y: 12462});
  data.push({t: moment('12-03-2020', date_format), y: 15113});
  data.push({t: moment('13-03-2020', date_format), y: 17660});
  data.push({t: moment('14-03-2020', date_format), y: 21157});
  data.push({t: moment('15-03-2020', date_format), y: 24747});
  data.push({t: moment('16-03-2020', date_format), y: 27980});
  data.push({t: moment('17-03-2020', date_format), y: 31506});
  data.push({t: moment('18-03-2020', date_format), y: 35713});
  data.push({t: moment('19-03-2020', date_format), y: 41035});
  data.push({t: moment('20-03-2020', date_format), y: 47021});
  data.push({t: moment('21-03-2020', date_format), y: 53578});
  data.push({t: moment('22-03-2020', date_format), y: 59138});
  return data;
}

var data_raw_SL = getRawDataSriLanka();
var data_SL = getDataSriLanka();
var data_IT = getDataItaly();

var sim_params = initializeSimulationParameters(data_SL.length);

var data_predicted = getPredictionData(data_SL[0].t);

var chart = [];
var chart_config = [];

window.onload = function()
{
  updateChart();

  let slider_interv0_T = document.getElementById('slider_interv0_T');
  let slider_interv1_T = document.getElementById('slider_interv1_T');
  slider_interv0_T.max = sim_params.T_hist + sim_params.T_pred;
  slider_interv1_T.max = sim_params.T_hist + sim_params.T_pred;
  slider_interv0_T.value = 18; //slider_interv1_T.max;
  slider_interv1_T.value = slider_interv0_T.max;
  document.getElementById('slider_interv0_T_value').innerHTML = slider_interv0_T.value;
  document.getElementById('slider_interv1_T_value').innerHTML = slider_interv1_T.value;
}

function formatNumber(num) {
  return num.toString().replace(/(\d)(?=(\d{3})+(?!\d))/g, '$1,')
}

function updateChart()
{
  let canvas = document.getElementById('chart_canvas');
  //canvas.width = 0.7*window.innerWidth;
  canvas.height = 0.65*window.innerHeight;
  let ctx = canvas.getContext('2d');

  let check_logy = document.getElementById('check_log_y');

  let xaxis_config = {
	    type: 'time',
			distribution: 'linear',
			time: {
			  tooltipFormat: 'MMM D',
			  unit: 'day'
//        displayFormats: {
//          day: 'DD-MM-YYYY'
//        }
      },
			offset: true,
			scaleLabel: {
			  display: true,
			  labelString: 'Date',
			  fontSize: 15
			},
			stacked: true
		};

  let yaxis_config = {
			gridLines: {
				drawBorder: false
			},
			type: 'linear',
			scaleLabel: {
				display: true,
				labelString: 'No. of cases',
				fontSize: 15
			}
		};

  if (check_logy.checked)
  {
    yaxis_config.type = 'logarithmic';
    yaxis_config.ticks = logarithmic_ticks;
  }

  const bar_width_frac = 1.0;
  const cat_width_frac = 0.9;

	chart_config = {
			data: {
				datasets: [
				{
					label: 'Susceptible',
					backgroundColor: '#eeeeec',
					data: data_predicted.categorized[0],
					type: 'bar',
					stack: 'stack0',
					hidden: true,
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 2
				},
				{
					label: 'Exposed',
					backgroundColor: 'rgba(255, 220, 160, 0.75)', //'#e9b96e',
					data: data_predicted.categorized[1],
					type: 'bar',
					stack: 'stack0',
					hidden: true,
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 3
				},
				{
					label: 'Mildly infected',
					backgroundColor: 'rgba(255, 200, 0, 0.75)',
					data: data_predicted.categorized[2],
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 4
				},
				{
					label: 'Severely infected',
					backgroundColor: 'rgba(240, 150, 40, 0.75)',
					data: data_predicted.categorized[3],
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 5
				},
				{
					label: 'Critically infected',
					backgroundColor: 'rgba(200, 0, 0, 0.75)',
					data: data_predicted.categorized[4],
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 6
				},
				{
					label: 'Recovered',
					backgroundColor: 'rgba(130, 210, 50, 0.75)', //#73d216',
					data: data_predicted.categorized[5],
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 7
				},
				{
					label: 'Fatal',
					backgroundColor: 'rgba(10, 10, 10, 0.75)',
					data: data_predicted.categorized[6],
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 8
				},
				{
					label: 'Sri Lanka - actual diagnosed',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#3465a4',
					data: data_SL,
					type: 'line',
					fill: true,
					borderWidth: 2,
					order: 0
				},
				{
					label: 'Sri Lanka - predicted diagnosed',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#729fcf',
					borderDash: [5, 5],
					data: data_predicted.aggregated,
					type: 'line',
					fill: true,
					borderWidth: 2,
					order: 1
				},
//				{
//					label: 'Italy - confirmed cases',
//					backgroundColor: 'rgba(1,1,1,0)',
//					borderColor: '#fcaf3e',
//					data: data_IT,
//					type: 'line',
//					fill: true,
//					borderWidth: 2,
//					hidden: true
//				}
				]
			},
			options: {
			  responsive: false,
        maintainAspectRatio: false,
        animation: {
            duration: 500 // general animation time
        },
				scales: {
					xAxes: [xaxis_config],
					yAxes: [yaxis_config]
				},
				legend: {
            display: true,
            boxWidth: 10,
        },
        tooltips: {
          callbacks: {
                label: function(tooltipItem, data) {
                    //var type = data.datasets[tooltipItem.datasetIndex].label;
                    //var value = data.datasets[tooltipItem.datasetIndex].data[tooltipItem.index].y;

                    if (tooltipItem.datasetIndex >= 7)
                      return (data.datasets[tooltipItem.datasetIndex].label +
                              ": " + formatNumber(data.datasets[tooltipItem.datasetIndex].data[tooltipItem.index].y));

                    // Loop through all datasets to get the actual total of the index
                    var total = 0;
                    let labels = [];
                    for (var i = 2; i < 7; i++)
                    {
                      labels.push(data.datasets[i].label + ": " + formatNumber(data.datasets[i].data[tooltipItem.index].y));
                      total += data.datasets[i].data[tooltipItem.index].y;
                    }
                    labels.push("Total: " + formatNumber(total));
                    return labels;
                }
            }
        },
        annotation: {
          drawTime: 'afterDatasetsDraw',
          // Array of annotation configuration objects
          // See below for detailed descriptions of the annotation options
          annotations: [{
            drawTime: 'afterDraw', // overrides annotation.drawTime if set
            id: 'vertline_T0', // optional
            type: 'line',
            mode: 'vertical',
            scaleID: 'x-axis-0',
            value: data_predicted.aggregated[18-1].t,
            borderColor: 'rgba(50,50,50,0.5)',
            borderWidth: 2,
            borderDash: [2, 2]
          },
          {
            drawTime: 'afterDraw', // overrides annotation.drawTime if set
            id: 'vertline_T1', // optional
            type: 'line',
            mode: 'vertical',
            scaleID: 'x-axis-0',
            value: data_predicted.aggregated[data_predicted.aggregated.length-1].t,
            borderColor: 'rgba(50,50,50,0.0)',
            borderWidth: 2,
            borderDash: [2, 2],
            hidden: true
          }]
        },
				plugins: {
	        zoom: {
		        // Container for pan options
		        pan: {
			        enabled: true,
			        mode: 'x',
			        rangeMin: { x: null, y: 0 },
			        rangeMax: { x: null, y: null },
			        speed: 20,		// On category scale, factor of pan velocity
			        threshold: 10, // Minimal pan distance required before actually applying pan
		        },

		        // Container for zoom options
		        zoom: {
			        enabled: true,
			        drag: false, // Enable drag-to-zoom behavior
			        mode: 'x',
			        rangeMin: { x: null, y: 0 },
			        rangeMax: { x: null, y: null },
			        speed: 0.1, // (percentage of zoom on a wheel event)
			        sensitivity: 3, // On category scale, minimal zoom level before actually applying zoom
		        }
	        }
        }

			}
		};

  chart = new Chart(ctx, chart_config);
};

function setLogYAxis(is_log)
{
  let logarithmic_ticks = {
    min: 0,
    //max: 100000,
    callback: function (value, index, values) {
      if (value === 10000000) return "10M";
      if (value === 1000000) return "1M";
      if (value === 100000) return "100k";
      if (value === 10000) return "10k";
      if (value === 1000) return "1k";
      if (value === 100) return "100";
      if (value === 10) return "10";
      if (value === 0) return "0";
      return null;
    }
  };

  if (is_log)
  {
    chart_config.options.scales.yAxes[0].type = 'logarithmic';
    chart_config.options.scales.yAxes[0].ticks = logarithmic_ticks;
  }
  else
  {
    chart_config.options.scales.yAxes[0].type = 'linear';
    chart_config.options.scales.yAxes[0].ticks = {display: true};
  }
  chart.update();
};

function setDataOverlayItaly(overlay)
{
  chart_config.data.datasets[9].hidden = !overlay;
  chart.update();
}

function alignTimelines(align)
{
  if (align)
  {
    for (let datapair of chart_config.data.datasets[9].data)
      datapair.t.add(41,'days'); //Add offset to Italy dataset
  }
  else //load original datasets
  {
    chart_config.data.datasets[9].data = getDataItaly();
  }
  chart.update();
}

function resetZoom()
{
  chart.resetZoom();
}

function initializeSimulationParameters(hist_length)
{
  const pred_length = 7; //no. of days to predict
  const total_length = hist_length + 200; //pred_length;

  let q_input = new Array(total_length).fill(0.0);
  for (let i = 0; i < data_raw_SL.length; ++i)
    q_input[i] = data_raw_SL[i].y[3];

  let params = {
    T_hist: hist_length,
    T_pred: pred_length,
    dt: 0.5/24.0,                           //timestep size [days]
    b1N: new Array(total_length).fill(0.6), //transmission rate from mild to susceptible
    b2N: new Array(total_length).fill(0.0), //transmission rate from severe to susceptible
    b3N: new Array(total_length).fill(0.0), //transmission rate from critical to susceptible
    quarantine_input: q_input,              //no. of patients added directly to quarantine
    diag_frac: 0.7,                         //fraction of I1 patients that are diagnosed
    E_0: 5,                                 //number of individiuals exposed at start
    I1h_0: 0                                //number of hidden mildly infected people at start
  }

  for (let i = 18; i < total_length; ++i)
    params.b1N[i] = 0.1;

  return params;
}

function updateInterventions(ind)
{
  let container = document.getElementById("container_interv" + ind);

  if (document.getElementById("check_interv" + ind).checked)
  {
    container.style.color = "black";
    document.getElementById("slider_interv" + ind + "_T").disabled = false;
    document.getElementById("slider_interv" + ind + "_b1").disabled = false;
    chart.annotation.elements["vertline_T" + ind].options.borderColor = "rgba(50,50,50,0.5)";
  }
  else
  {
    container.style.color = "grey";
    document.getElementById("slider_interv" + ind + "_T").disabled = true;
    document.getElementById("slider_interv" + ind + "_b1").disabled = true;
    chart.annotation.elements["vertline_T" + ind].options.borderColor = "rgba(50,50,50,0.0)";
  }
  updateParameters(true);
}

function updateParameters(force = false)
{
  let requires_update = force;

  let b1 = Number(document.getElementById("slider_b1").value);
  document.getElementById("slider_b1_value").innerHTML = b1.toFixed(2);

  let interv_params = [{t: Infinity, val: 0}, {t: Infinity, val: 0}];

  for (let i = 0; i < 2; ++i)
  {
    if (document.getElementById("check_interv" + i).checked)
    {
      interv_params[i].t = Number(document.getElementById("slider_interv" + i + "_T").value);
      interv_params[i].val = Number(document.getElementById("slider_interv" + i + "_b1").value);
      document.getElementById("slider_interv" + i + "_T_value").innerHTML = interv_params[i].t;
      document.getElementById("slider_interv" + i + "_b1_value").innerHTML = interv_params[i].val.toFixed(2);
    }
  }

  for (let j = 0; j < sim_params.b1N.length; ++j)
  {
    let val = b1;

    let largest_t = -1;
    for (let k = 0; k < 2; ++k)
    {
      if (j >= interv_params[k].t && interv_params[k].t > largest_t)
      {
        largest_t = interv_params[k].t;
        val = interv_params[k].val;
      }
    }

    if (sim_params.b1N[j] != val)
    {
      sim_params.b1N[j] = val;
      requires_update = true;
    }
  }

  let slider_element_ids = ["slider_b2", "slider_b3"];
  let param_arrays = [sim_params.b2N, sim_params.b3N];

  for (let i = 0; i < 2; ++i)
  {
    let slider = document.getElementById(slider_element_ids[i]);
    if (slider)
    {
      let val = Number(slider.value);
      for (let j = 0; j < param_arrays[i].length; ++j)
        if (param_arrays[i][j] != val)
        {
          param_arrays[i][j] = val;
          requires_update = true;
          document.getElementById(slider_element_ids[i] + "_value").innerHTML = val.toFixed(2);
        }
    }
  }

  let slider_finalT = document.getElementById("slider_finalT");
  if (slider_finalT)
  {
    let val = Number(slider_finalT.value);
    if (sim_params.T_pred != val)
    {
      sim_params.T_pred = val;
      requires_update = true;
      document.getElementById("slider_finalT_text").innerHTML = "Predict for " + val + " days";
      slider_interv0_T.max = sim_params.T_hist + sim_params.T_pred;
      slider_interv1_T.max = sim_params.T_hist + sim_params.T_pred;
    }
  }

  if (requires_update)
  {
    data_predicted = getPredictionData(data_SL[0].t);
    for (let i = 0; i < 7; ++i)
      chart_config.data.datasets[i].data = data_predicted.categorized[i];

    chart_config.data.datasets[8].data = data_predicted.aggregated;

    for (let i = 0; i < 2; ++i)
    {
      if (interv_params[i].t != Infinity)
        chart.annotation.elements["vertline_T" + i].options.value = data_predicted.aggregated[interv_params[i].t-1].t;
    }
    chart.update();
  }
}

function getPredictionData(start_date)
{
  let sol_history = predictModel(sim_params);

  let data_agg = [];
  let data_cat = new Array(7);
  for (let i = 0; i < 7; ++i)
    data_cat[i] = new Array();

  const report_sum_indices = [2, 4, 5, 6, 7, 9];

  for (let i = 0; i < sol_history.length; i++)
  {
    let date = start_date.clone().add(i,'days');

    //Accumulate data into categories for plotting
    data_cat[0].push({t: date, y: Math.round(sol_history[i][0])}); //S
    data_cat[1].push({t: date, y: Math.round(sol_history[i][1])}); //E
    data_cat[2].push({t: date, y: Math.round(sol_history[i][2] + sol_history[i][3] + sol_history[i][4])}); //I1
    data_cat[3].push({t: date, y: Math.round(sol_history[i][5])}); //I2
    data_cat[4].push({t: date, y: Math.round(sol_history[i][6])}); //I3
    data_cat[5].push({t: date, y: Math.round(sol_history[i][7] + sol_history[i][8])}); //R
    data_cat[6].push({t: date, y: Math.round(sol_history[i][9])}); //D

    let num_confirmed_cases = 0;
    for (let j = 0; j < report_sum_indices.length; ++j)
      num_confirmed_cases += sol_history[i][report_sum_indices[j]];
    data_agg.push({t: date, y: Math.round(num_confirmed_cases)});
  }

  return {aggregated: data_agg, categorized: data_cat};
}

function predictModel(params)
{
  //Periods [days]
  let T_incub  = 5;
  let T_mild   = 6;
  let T_severe = 4;
  let T_icu    = 10;

  //Probabilities
  let prob_I1_E   = 1;  //exposed to mildly infected
  let prob_R_I1   = 0.81*prob_I1_E; //mild to recovered
  let prob_I2_I1  = 1 - prob_R_I1; //mild to severe
  let prob_R_I2   = 0.14/prob_I2_I1; //severe to recovered
  let prob_I3_I2  = 1 - prob_R_I2; //severe to critical
  let prob_D_I3   = 0.02/(prob_I3_I2*prob_I2_I1); //critical to dead
  let prob_R_I3   = 1 - prob_D_I3; //critical to recovered

  //Rates [1/day]
  let a   = (1/T_incub)  * prob_I1_E;
  let g1  = (1/T_mild)   * prob_R_I1;
  let p1  = (1/T_mild)   * prob_I2_I1;
  let g2  = (1/T_severe) * prob_R_I2;
  let p2  = (1/T_severe) * prob_I3_I2;
  let g3  = (1/T_icu)    * prob_R_I3;
  let mu  = (1/T_icu)    * prob_D_I3;

  let N = 21.4e6; //Population of Sri Lanka

  //Initial solution
  let E_0 = params.E_0;
  let I1d_0 = 0;
  let I1h_0 = params.I1h_0;
  let I1q_0 = 0;
  let I2_0 = 0;
  let I3_0 = 0;
  let Rd_0 = 1;  //One patient had already recovered in SL
  let Rh_0 = 0;
  let D_0 = 0;
  let S_0 = N - E_0 - I1d_0 - I1h_0 - I1q_0 - I2_0 - I3_0 - Rd_0 - Rh_0 - D_0;

  //Solution vector: [S, E, I1d, I1h, I1q, I2, I3, Rd, Rh, D]
  let solution_hist = [[S_0, E_0, I1d_0, I1h_0, I1q_0, I2_0, I3_0, Rd_0, Rh_0, D_0]];

  const nt = params.T_hist + params.T_pred - 1;
  const nt_sub = 1.0/params.dt;
  const c = params.diag_frac;

  for (let i = 0; i < nt; i++)
  {
    let sol = solution_hist[i].slice(); //copy last solution vector

    let b1 = params.b1N[i] / N;
    let b2 = params.b2N[i] / N;
    let b3 = params.b3N[i] / N;
    let q_input = params.quarantine_input[i];

    for (let j = 0; j < nt_sub; j++)
    {
      let dS = -(b1*(sol[2] + sol[3]) + b2*sol[5] + b3*sol[6])*sol[0];
      let dsol = [dS,                                               //S
                  -dS - a*sol[1],                                   //E
                  a*c*sol[1]     - (g1 + p1)*sol[2],                //I1d
                  a*(1-c)*sol[1] - (g1 + p1)*sol[3],                //I1h
                  q_input        - (g1 + p1)*sol[4],                //I1q
                  p1*(sol[2] + sol[3] + sol[4]) - (g2 + p2)*sol[5], //I2
                  p2*sol[5] - (g3 + mu)*sol[6],                     //I3
                  g1*(sol[2] + sol[4]) + g2*sol[5] + g3*sol[6],     //Rd
                  g1*(sol[3]),                                      //Rh
                  mu*sol[6]                                         //D
                 ];

      for (let k = 0; k < sol.length; k++)
        sol[k] += dsol[k]*params.dt;
    } //sub-timestepping [hrs]

    if (sol[1] < 0.5)
      sol[1] = 0.0;

    solution_hist.push(sol); //save solution daily
  }

  return solution_hist;
}


function optimizeParameters()
{
  let params = initializeSimulationParameters(data_SL.length);
  let T_pred_init = params.T_pred;
  params.T_pred = 0;

  let res = getFitResidual(params);
  let resnorm = getL2Norm(res);

  let m = res.length;
  let n_b1 = (params.T_hist - 1); //no. of beta_1 values to optimize
  let n = n_b1 + 2; //beta_1 values, E_0, I1h_0

  let param_vec = new Array(n).fill(0);
  for (let i = 0; i < n_b1; ++i)
    param_vec[i] = params.b1N[i];
  param_vec[n_b1] = params.E_0;
  param_vec[n_b1 + 1] = params.I1h_0;

  let dparam_vec = new Array(n).fill(0);
  let param_vec0 = new Array(n).fill(0);

  for (let iter = 0; iter < 100; ++iter)
  {
    console.log("Iter " + iter + ": " + resnorm);

    let jac = getFitJacobian(params);

    //Update parameter vector using gradient descent:
    //u(n+1) = u(n) - eta * (dR/du)^T R(u(n))
    for (let i = 0; i < n; ++i)
    {
      let dp = 0.0;
      for (let j = 0; j < m; ++j)
        dp += jac[j*n + i] * res[j]
      dparam_vec[i] = dp;
      param_vec0[i] = param_vec[i];
    }

    let eta = 0.5;

    while (eta >= 1e-7)
    {
      for (let i = 0; i < n; ++i)
      {
        param_vec[i] = param_vec0[i] - eta*dparam_vec[i];

        if (i < n_b1)
        {
          if (param_vec[i] < 0.0)
            param_vec[i] = 0.0;
          params.b1N[i] = param_vec[i];
        }
        else if (i == n_b1)
        {
          if (param_vec[i] < 0.0)
            param_vec[i] = 0.0;
          params.E_0 = param_vec[i];
        }
        else if (i == n_b1+1)
        {
          if (param_vec[i] < 0.0)
            param_vec[i] = 0.0;
          params.I1h_0 = param_vec[i];
        }
      }

      //Evaluate new residual
      res = getFitResidual(params);
      let resnorm_new = getL2Norm(res);
      //console.log("  " + eta + ", " + getL2Norm(res));

      if (resnorm_new < resnorm)
      {
        resnorm = resnorm_new;
        break;
      }

      eta /= 2.0;
    } //linesearch
  } //gradient descent

  for (let i = 0; i < n_b1; ++i)
    sim_params.b1N[i] = params.b1N[i];
  sim_params.E_0 = params.E_0;
  sim_params.I1h_0 = params.I1h_0;

  for (let i = n_b1; i < sim_params.b1N.length; ++i)
    sim_params.b1N[i] = sim_params.b1N[n_b1-1];

  updateParameters();
  return params;
}

function getFitResidual(params)
{
  const num_eq = 2;
  let sol_hist = predictModel(params);
  let residual = new Array(num_eq*sol_hist.length).fill(0);

  for (let i = 0; i < sol_hist.length; ++i)
  {
    //I1d + I1q + I2 + I3 = no. of diagnosed patients in hospitals
    residual[num_eq*i] = sol_hist[i][2] + sol_hist[i][4] + sol_hist[i][5] + sol_hist[i][6] - data_raw_SL[i].y[0];
    //residual[num_eq*i] = sol_hist[i][2] + sol_hist[i][4] + sol_hist[i][5] + sol_hist[i][6] + sol_hist[i][7]
    //                     - data_raw_SL[i].y[0] - data_raw_SL[i].y[1];

    //Rd = no. of recovered patients
    residual[num_eq*i + 1] = sol_hist[i][7] - data_raw_SL[i].y[1];

    //D = no. of fatalities
    //residual[num_eq*i + 2] = sol_hist[i][9] - data_raw_SL[i].y[2];
  }
  return residual;
}

function getFitJacobian(params)
{
  let m = 2*params.T_hist;
  let n_b1 = (params.T_hist - 1); //no. of beta_1 values to optimize
  let n = n_b1 + 2; //beta_1 values, E_0, I1h_0
  let jac = new Array(m*n).fill(0);

  const delta_b1 = 1e-5;
  const delta_E0 = 1e-5;
  const delta_I1h0 = 1e-5;

  for (let j = 0; j < n_b1; ++j)
  {
    //Compute finite difference
    params.b1N[j] += delta_b1;
    let Rp = getFitResidual(params);
    params.b1N[j] -= 2*delta_b1;
    let Rm = getFitResidual(params);
    params.b1N[j] += delta_b1;

    for (let i = 0; i < m; ++i)
      jac[i*n + j] = (Rp[i] - Rm[i])/(2*delta_b1);
  }

  params.E_0 += delta_E0;
  let Rp = getFitResidual(params);
  params.E_0 -= 2*delta_E0;
  let Rm = getFitResidual(params);
  params.E_0 += delta_E0;

  for (let i = 0; i < m; ++i)
    jac[i*n + n_b1] = (Rp[i] - Rm[i])/(2*delta_E0);

  params.I1h_0 += delta_I1h0;
  Rp = getFitResidual(params);
  params.I1h_0 -= 2*delta_I1h0;
  Rm = getFitResidual(params);
  params.I1h_0 += delta_I1h0;

  for (let i = 0; i < m; ++i)
    jac[i*n + n_b1 + 1] = (Rp[i] - Rm[i])/(2*delta_I1h0);

  return jac;
}

function getL2Norm(vec)
{
  let norm = 0.0;
  for (let i = 0; i < vec.length; ++i)
    norm += vec[i]*vec[i];
  return Math.sqrt(norm);
}
