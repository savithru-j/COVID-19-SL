// COVID-19 Simulation Tool in JavaScript
// Copyright 2020 Savithru Jayasinghe and Dushan Wadduwage
// Licensed under the MIT License (LICENSE.txt)

'use strict';

// var date_format = 'DD-MM-YYYY';
var date_format = 'YYYY-MM-DD';

var population_data = {
  "Sri Lanka" : 21.44E+6,
  "Singapore" : 5.61E+6,
  "US" : 327.2E+6,
  "United Kingdom" : 66.65E+6,
}

//Data about any quarantined index individuals in each country (required for Iq category)
var quarantine_data = {
  "Sri Lanka" : [{t: "2020-03-13", y: 2},
                 {t: "2020-03-14", y: 2},
                 {t: "2020-03-15", y: 7},
                 {t: "2020-03-16", y: 6},
                 {t: "2020-03-17", y: 4},
                 {t: "2020-03-18", y: 1},
                 {t: "2020-03-19", y: 2},
                 {t: "2020-03-20", y: 3},
                 {t: "2020-03-21", y: 6}]
}

var data_start_dates = {
    "Sri Lanka": "2020-03-01"
};

var custom_country_data = {
  "Sri Lanka" : {
      t_start: [0, 13], //indices to start dates of any interventions
      b1N: [0.85, 0.125], //values of b1N for each intervention segment defined in t_start
      b2N: [0, 0], //values of b2N
      b3N: [0, 0],
      diag_frac: [0.11, 0.11],
      E0_0: 5, //no. of individuals exposed at start
      Rd_0: 1, //no. of recovered-diagnosed individuals at start
  }
}

//The control parameters will be set to these default values when the user first loads the page.
var default_controls = {
  T_pred: 7,    //prediction length
  b1N: 0.0,     //beta_1 value
  b2N: 0.0,     //beta_2 value
  b3N: 0.0,     //beta_3 value
  diag_frac: 0.5 //fraction of mild-patients that are diagnosed
}

//Storage for real and predicted data to be plotted on chart: {categorized, total}
var data_real;
var data_predicted;
var sim_params;

//Chart objects
var main_chart, control_chart;
var control_chart_active_point = null;
var control_chart_canvas = null;
var last_active_tooltip_day = 0;
var active_control_parameter = "b1N";
var active_country = "Sri Lanka";


window.onload = function()
{
  generateCountryDropDown();
  changeCountry(active_country);

  //Update UI controls to match default values
  document.getElementById("slider_finalT").value = default_controls.T_pred;
  document.getElementById("slider_finalT_value").innerHTML = default_controls.T_pred;
}

function generateCountryDropDown()
{
  let menu = document.getElementById("dropdown_country");
  let country_list = [];
  for (let key of Object.keys(world_data))
    country_list.push(key);
  country_list.sort();

  for (let name of country_list)
  {
    let option = document.createElement("option");
    option.text = name;
    menu.add(option);
  }
  menu.value = active_country;
}

function changeCountry(country_name)
{
  active_country = country_name;
  data_real = getCountryData(country_name);
  sim_params = initializeSimulationParameters(data_real.total.length, default_controls.T_pred);

  data_predicted = getPredictionData(data_real.total[0].t);
  document.getElementById("prediction_error").innerHTML = getCurrentPredictionError().toFixed(3);

  if (main_chart)
    refreshAllChartData();
  else {
    setupChart();
    updateLegend();
    setupControlChart();
  }
}

function getCountryData(country_name)
{
  let data_array = world_data[country_name];

  let start_date = data_start_dates[country_name];
  if (start_date)
    start_date = moment(start_date, date_format);

  let data_cat = [], data_total = [];
  if (start_date)
  {
    for (let data of data_array)
    {
      let data_t = moment(data.date, date_format);
      if (start_date <= data_t)
      {
        data_total.push({t: data_t, y: data.confirmed});
        data_cat.push({t: data_t, y: [data.confirmed, data.recovered, data.deaths, 0]});
      }
    }
  }
  else
  {
    for (let data of data_array)
    {
      let data_t = moment(data.date, date_format);
      if (data.confirmed + data.recovered + data.deaths > 0)
      {
        data_total.push({t: data_t, y: data.confirmed});
        data_cat.push({t: data_t, y: [data.confirmed, data.recovered, data.deaths, 0]});
      }
    }
  }

  //Check if country has any quarantine data
  let qdata = quarantine_data[country_name];
  if (qdata)
  {
    for (let data of qdata)
    {
      let data_t = moment(data.t, date_format);
      let ind = data_t.diff(data_cat[0].t, 'days');
      data_cat[ind].y[3] = data.y;
    }
  }

  return {total: data_total, categorized: data_cat};
}

function customizeParametersByCountry(country_name, params)
{
  let data = custom_country_data[country_name];
  if (data)
  {
    for (let i = 0; i < data.t_start.length; ++i)
    {
      let j_end = (i < data.t_start.length-1) ? data.t_start[i+1] : params.b1N.length;
      for (let j = data.t_start[i]; j < j_end; ++j)
      {
        params.b1N[j] = data.b1N[i];
        params.b2N[j] = data.b2N[i];
        params.b3N[j] = data.b3N[i];
        params.diag_frac[j] = data.diag_frac[i];
      }
    }

    //Update initial E0 and Ru
    if (data.E0_0)
      params.E0_0 = data.E0_0;

    if (data.Rd_0)
      params.Rd_0 = data.Rd_0;
  }

  //Update quarantine input data
  for (let i = 0; i < data_real.categorized.length; ++i)
    params.quarantine_input[i] = data_real.categorized[i].y[3];

  //Update population
  let N = population_data[country_name];
  if (N)
    params.population = N;
  else
    console.log("Population data not found!");
}

function formatNumber(num) {
  return num.toString().replace(/(\d)(?=(\d{3})+(?!\d))/g, '$1,')
}

function setupChart()
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

	let main_chart_config = {
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
					order: 11
				},
				{
					label: 'Exposed',
					backgroundColor: 'rgba(255, 220, 160, 0.75)',
					data: data_predicted.categorized[1],
					type: 'bar',
					stack: 'stack0',
					hidden: true,
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 10
				},
        {
          label: 'Asymptomatic',
          backgroundColor: 'rgba(180, 240, 255, 0.75)',
          data: data_predicted.categorized[2],
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 9
        },
        {
          label: 'Mildly infected - unreported',
          backgroundColor: 'rgba(255, 200, 0, 0.5)',
          data: data_predicted.categorized[3],
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 8
        },
        {
          label: 'Recovered - unreported',
          backgroundColor: 'rgba(130, 210, 50, 0.4)',
          data: data_predicted.categorized[4],
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 7
        },
        {
          label: 'Recovered',
          backgroundColor: 'rgba(130, 210, 50, 0.75)', //#73d216',
          data: data_predicted.categorized[5],
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 6
        },
        {
          label: 'Fatal',
          backgroundColor: 'rgba(10, 10, 10, 0.75)',
          data: data_predicted.categorized[6],
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 5
        },
        {
          label: 'Critically infected',
          backgroundColor: 'rgba(200, 0, 0, 0.75)',
          data: data_predicted.categorized[7],
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 4
        },
        {
          label: 'Severely infected',
          backgroundColor: 'rgba(240, 150, 40, 0.75)',
          data: data_predicted.categorized[8],
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 3
        },
				{
					label: 'Mildly infected',
					backgroundColor: 'rgba(255, 200, 0, 0.75)',
					data: data_predicted.categorized[9],
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 2
				},
				{
					label: 'Actual diagnosed',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#3465a4',
					data: data_real.total,
					type: 'line',
					fill: true,
					borderWidth: 2,
					order: 0
				},
				{
					label: 'Predicted diagnosed',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#729fcf',
					borderDash: [5, 5],
					data: data_predicted.total,
					type: 'line',
					fill: true,
					borderWidth: 2,
					order: 1
				}
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
            boxWidth: 10
        },
        tooltips: {
          callbacks: {
                label: function(tooltipItem, data) {
                    //var type = data.datasets[tooltipItem.datasetIndex].label;
                    //var value = data.datasets[tooltipItem.datasetIndex].data[tooltipItem.index].y;

                    updateLegend(tooltipItem.index);

                    if (tooltipItem.datasetIndex >= 10)
                      return (data.datasets[tooltipItem.datasetIndex].label +
                              ": " + formatNumber(data.datasets[tooltipItem.datasetIndex].data[tooltipItem.index].y));

                    // Loop through all datasets to get the actual total of the index
                    // var total = 0;
                    // let labels = [];
                    // for (var i = 2; i < 10; i++)
                    // {
                    //   labels.push(data.datasets[i].label + ": " + formatNumber(data.datasets[i].data[tooltipItem.index].y));
                    //   total += data.datasets[i].data[tooltipItem.index].y;
                    // }
                    // labels.push("Total: " + formatNumber(total));
                    // return labels;
                }
            }
        },
        // annotation: {
        //   drawTime: 'afterDatasetsDraw',
        //   // Array of annotation configuration objects
        //   // See below for detailed descriptions of the annotation options
        //   annotations: [{
        //     drawTime: 'afterDraw', // overrides annotation.drawTime if set
        //     id: 'vertline_T0', // optional
        //     type: 'line',
        //     mode: 'vertical',
        //     scaleID: 'x-axis-0',
        //     value: data_predicted.total[default_controls.T0-1].t,
        //     borderColor: 'rgba(50,50,50,0.5)',
        //     borderWidth: 2,
        //     borderDash: [2, 2]
        //   },
        //   {
        //     drawTime: 'afterDraw', // overrides annotation.drawTime if set
        //     id: 'vertline_T1', // optional
        //     type: 'line',
        //     mode: 'vertical',
        //     scaleID: 'x-axis-0',
        //     value: data_predicted.total[data_predicted.total.length-1].t,
        //     borderColor: 'rgba(50,50,50,0.0)',
        //     borderWidth: 2,
        //     borderDash: [2, 2],
        //     hidden: true
        //   }]
        // },
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
              onPan: function () { syncPanAndZoom(main_chart, control_chart); }
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
              onZoom: function () { syncPanAndZoom(main_chart, control_chart); }
		        }
	        }
        }

			}
		};

  main_chart = new Chart(ctx, main_chart_config);

  //Save off original axes (before any pan/zoom is applied)
  main_chart.$zoom._originalOptions[main_chart.options.scales.xAxes[0].id] = main_chart.options.scales.xAxes[0];
  main_chart.$zoom._originalOptions[main_chart.options.scales.yAxes[0].id] = main_chart.options.scales.yAxes[0];
};

function setupControlChart()
{
  control_chart_canvas = document.getElementById('control_chart_canvas');
  control_chart_canvas.height = 0.2*window.innerHeight;
  let ctx = control_chart_canvas.getContext('2d');

  let xaxis_config = {
      type: 'time',
      distribution: 'linear',
      time: {
        tooltipFormat: 'MMM D',
        unit: 'day'
      },
      offset: true,
      scaleLabel: {
        display: false
      },
      ticks: {
          display: false //this will remove only the label
      }
    };

  let yaxis_config = {
      gridLines: {
        drawBorder: false
      },
      type: 'linear',
      scaleLabel: {
        display: true,
        labelString: 'Beta_1',
        fontSize: 15
      },
      ticks: {
          min: 0,
          max: 1
      }
    };

  let control_chart_config = {
      data: {
        datasets: [{
          label: 'beta_1',
          backgroundColor: 'rgba(1,1,1,0)',
          borderColor: 'rgb(50, 160, 220)',
          data: getControlChartData(),
          type: 'line',
          fill: true,
          borderWidth: 2,
          order: 1,
          cubicInterpolationMode: 'monotone'
        }]
      },
      options: {
        responsive: false,
        maintainAspectRatio: false,
        animation: {
            duration: 200 // general animation time
        },
        scales: {
          xAxes: [xaxis_config],
          yAxes: [yaxis_config]
        },
        legend: {
            display: false,
        },
        hover: {
          onHover: function(e, el) {
            control_chart_canvas.style.cursor = el[0] ? "pointer" : "default";
          }
        },
        plugins: {
          zoom: {
            // Container for pan options
            pan: {
              enabled: true,
              mode: 'x',
              rangeMin: { x: null, y: 0 },
              rangeMax: { x: null, y: 1 },
              speed: 20,		// On category scale, factor of pan velocity
              threshold: 10, // Minimal pan distance required before actually applying pan
              onPan: function () { syncPanAndZoom(control_chart, main_chart); }
            },

            // Container for zoom options
            zoom: {
              enabled: true,
              drag: false, // Enable drag-to-zoom behavior
              mode: 'x',
              rangeMin: { x: null, y: 0 },
              rangeMax: { x: null, y: 1 },
              speed: 0.1, // (percentage of zoom on a wheel event)
              sensitivity: 3, // On category scale, minimal zoom level before actually applying zoom
              onZoom: function () { syncPanAndZoom(control_chart, main_chart); }
            },
          }
        }

      }
    };
  control_chart = new Chart(ctx, control_chart_config);

  // set pointer event handlers for canvas element
  control_chart_canvas.onpointerdown = down_handler;
  control_chart_canvas.onpointerup = up_handler;
  control_chart_canvas.onpointermove = null;

  //Save off original axes (before any pan/zoom is applied)
  control_chart.$zoom._originalOptions[control_chart.options.scales.xAxes[0].id] = control_chart.options.scales.xAxes[0];
  control_chart.$zoom._originalOptions[control_chart.options.scales.yAxes[0].id] = control_chart.options.scales.yAxes[0];
}

function down_handler(event)
{
  // check for data point near event location
  const points = control_chart.getElementAtEvent(event, {intersect: false});
  if (points.length > 0) {
      // grab nearest point, start dragging
      control_chart_active_point = points[0];
      control_chart_canvas.onpointermove = move_handler;
  };
};

function up_handler(event)
{
  // release grabbed point, stop dragging
  control_chart_active_point = null;
  control_chart_canvas.onpointermove = null;
};

function move_handler(event)
{
  // locate grabbed point in chart data
  if (control_chart_active_point != null) {
      let data = control_chart_active_point._chart.data;
      let datasetIndex = control_chart_active_point._datasetIndex;

      // read mouse position
      const helpers = Chart.helpers;
      let position = helpers.getRelativePosition(event, control_chart);

      // convert mouse position to chart y axis value
      let chart_area = control_chart.chartArea;
      let yaxis = control_chart.scales["y-axis-0"];
      let yval_new = map(position.y, chart_area.bottom, chart_area.top, yaxis.min, yaxis.max);
      yval_new = Math.round(yval_new*1000)/1000;
      yval_new = Math.min(Math.max(yval_new, 0.0), 1.0);

      //Update values to the right of the current index, until a different value is encountered.
      let datavec = sim_params[active_control_parameter]; //get the data vector of the parameter current selected for editing
      let yval_old = datavec[control_chart_active_point._index];

      if (yval_new != yval_old)
      {
        for (let i = control_chart_active_point._index; i < datavec.length; ++i)
        {
          if (datavec[i] != yval_old)
            break;
          datavec[i] = yval_new;
        }

        updateParameters(true);
      }
  };
};

// map value to other coordinate system
function map(value, start1, stop1, start2, stop2) {
    return start2 + (stop2 - start2) * ((value - start1) / (stop1 - start1))
};

function syncPanAndZoom(chart_from, chart_to)
{
  chart_to.options.scales.xAxes[0].time.min = chart_from.options.scales.xAxes[0].time.min;
  chart_to.options.scales.xAxes[0].time.max = chart_from.options.scales.xAxes[0].time.max;
  chart_to.update();
};

function changeControlChartParameter(parameter)
{
  active_control_parameter = parameter;
  refreshControlChartData();
}

function getControlChartData()
{
  let data = [];
  for (let i = 0; i < data_predicted.total.length; ++i)
    data.push({t: data_predicted.total[i].t, y: sim_params[active_control_parameter][i]});
  return data;
}

function refreshAllChartData()
{
  refreshMainChartData();
  refreshControlChartData();
}

function refreshMainChartData()
{
  if (main_chart)
  {
    let n_cat = data_predicted.categorized.length;
    for (let i = 0; i < n_cat; ++i)
      main_chart.data.datasets[i].data = data_predicted.categorized[i];

    main_chart.data.datasets[n_cat].data = data_real.total;
    main_chart.data.datasets[n_cat+1].data = data_predicted.total;

    main_chart.update();
    delete main_chart.$zoom._originalOptions[main_chart.options.scales.xAxes[0].id].time.min;
    delete main_chart.$zoom._originalOptions[main_chart.options.scales.xAxes[0].id].time.max;

    updateLegend();
  }
}

function refreshControlChartData()
{
  const param_to_labelstring = {
    "b1N": "Beta_1",
    "b2N": "Beta_2",
    "b3N": "Beta_3",
    "diag_frac": "c"
  };

  if (control_chart)
  {
    control_chart.data.datasets[0].data = getControlChartData();
    control_chart.data.datasets[0].label = param_to_labelstring[active_control_parameter];
    control_chart.options.scales.yAxes[0].scaleLabel.labelString = param_to_labelstring[active_control_parameter];
    control_chart.update();
    delete control_chart.$zoom._originalOptions[control_chart.options.scales.xAxes[0].id].time.min;
    delete control_chart.$zoom._originalOptions[control_chart.options.scales.xAxes[0].id].time.max;
  }
}

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
    main_chart.options.scales.yAxes[0].type = 'logarithmic';
    main_chart.options.scales.yAxes[0].ticks = logarithmic_ticks;
  }
  else
  {
    main_chart.options.scales.yAxes[0].type = 'linear';
    main_chart.options.scales.yAxes[0].ticks = {display: true};
  }
  main_chart.update();
}

function showOnlyDiagnosed(flag)
{
  for (let i = 0; i < 5; ++i)
    main_chart.data.datasets[i].hidden = flag;
  if (!flag)
  { //Keep susceptible and exposed hidden
    main_chart.data.datasets[0].hidden = true;
    main_chart.data.datasets[1].hidden = true;
  }
  main_chart.update();
}

function resetZoom()
{
  main_chart.resetZoom();
  control_chart.resetZoom();
}

function updateLegend(day = last_active_tooltip_day)
{
  if (day < 0 || day >= data_predicted.total.length)
    return;

  last_active_tooltip_day = day;

  document.getElementById("legend_date").innerHTML = data_predicted.total[day].t.format("MMM-DD-YYYY");

  for (let i = 1; i < data_predicted.categorized.length; ++i)
    document.getElementById("legend_pred" + i).innerHTML = formatNumber(data_predicted.categorized[i][day].y);
  document.getElementById("legend_pred_infected").innerHTML = formatNumber(data_predicted.categorized[7][day].y + data_predicted.categorized[8][day].y + data_predicted.categorized[9][day].y);
  document.getElementById("legend_pred_total").innerHTML = formatNumber(data_predicted.total[day].y);

  let true_data = ['-', '-', '-', '-'];
  if (day < data_real.categorized.length)
  {
    true_data[0] = formatNumber(data_real.categorized[day].y[0] - data_real.categorized[day].y[1] - data_real.categorized[day].y[2]);
    true_data[1] = formatNumber(data_real.categorized[day].y[1]);
    true_data[2] = formatNumber(data_real.categorized[day].y[2]);
    true_data[3] = formatNumber(data_real.categorized[day].y[0]);
  }
  document.getElementById("legend_true_infected").innerHTML = true_data[0];
  document.getElementById("legend_true_recovered").innerHTML = true_data[1];
  document.getElementById("legend_true_fatal").innerHTML = true_data[2];
  document.getElementById("legend_true_total").innerHTML = true_data[3];
}

function initializeSimulationParameters(hist_length, pred_length)
{
  //allocate maximum possible through sliders so that we don't have to resize later
  let total_length = hist_length + 200;

  //Periods [days]
  const T_incub0 = 3;
  const T_incub1 = 2;
  const T_asympt = 6;
  const T_mild   = 6;
  const T_severe = 4;
  const T_icu    = 10;

  //Probabilities
  const prob_E0_E1 = 1;  //non-infectious exposed to infectious exposed
  const prob_E1_I0 = 0.3 / prob_E0_E1;  //exposed to asymptomatic
  const prob_E1_I1 = 1 - prob_E1_I0;  //exposed to mild
  const prob_I0_R  = 1;
  const prob_I1_R  = 0.8; //mild to recovered
  const prob_I1_I2 = 1 - prob_I1_R; //mild to severe
  const prob_I2_R  = 0.15/(prob_I1_I2); //severe to recovered
  const prob_I2_I3 = 1 - prob_I2_R; //severe to critical
  const prob_I3_D  = 0.02/(prob_I1_I2*prob_I2_I3); //critical to dead
  const prob_I3_R  = 1 - prob_I3_D; //critical to recovered

  let params = {
    T_hist: hist_length,
    T_pred: pred_length,
    dt: 0.5/24.0,                                                       //timestep size [days]
    b1N: new Array(total_length).fill(default_controls.b1N),            //transmission rate from mild to susceptible
    b2N: new Array(total_length).fill(default_controls.b2N),            //transmission rate from severe to susceptible
    b3N: new Array(total_length).fill(default_controls.b3N),            //transmission rate from critical to susceptible
    quarantine_input: new Array(total_length).fill(0.0),                //no. of patients added directly to quarantine
    diag_frac: new Array(total_length).fill(default_controls.diag_frac),//fraction of I1 patients that are diagnosed
    population: 1E7,                                                    //population of country
    E0_0: 5,                                                            //number of non-infectious exposed individuals at start
    Rd_0: 0,                                                            //number of recovered-diagnosed individuals at start
    //rate parameters below [1/day]
    a0: (1/T_incub0) * prob_E0_E1,
    a10: (1/T_incub1) * prob_E1_I0,
    a11: (1/T_incub1) * prob_E1_I1,
    g0: (1/T_asympt) * prob_I0_R,
    g1: (1/T_mild)   * prob_I1_R,
    p1: (1/T_mild)   * prob_I1_I2,
    g2: (1/T_severe) * prob_I2_R,
    p2: (1/T_severe) * prob_I2_I3,
    g3: (1/T_icu)    * prob_I3_R,
    mu: (1/T_icu)    * prob_I3_D,
  }

  //Modify any parameters that are specific to the currently active country.
  customizeParametersByCountry(active_country, params);

  return params;
}

function updateParameters(force = false)
{
  let requires_update = force;

  let slider_finalT = document.getElementById("slider_finalT");
  if (slider_finalT)
  {
    let val = Number(slider_finalT.value);
    if (sim_params.T_pred != val)
    {
      sim_params.T_pred = val;
      requires_update = true;
      document.getElementById("slider_finalT_value").innerHTML = val;
    }
  }

  if (requires_update)
  {
    data_predicted = getPredictionData(data_real.total[0].t);
    refreshAllChartData();
    document.getElementById("prediction_error").innerHTML = getCurrentPredictionError().toFixed(3);
  }
}

function getPredictionData(start_date)
{
  let sol_history = predictModel(sim_params);

  let data_agg = [];
  let data_cat = new Array(10);
  for (let i = 0; i < data_cat.length; ++i)
    data_cat[i] = new Array();

  const report_sum_indices = [4, 6, 7, 8, 9, 11]; //I1d + I1q + I2 + I3 + Rd + D

  for (let i = 0; i < sol_history.length; i++)
  {
    let date = start_date.clone().add(i,'days');
    let c = sim_params.diag_frac[i];

    //Accumulate data into categories for plotting
    data_cat[0].push({t: date, y: Math.round(sol_history[i][0])}); //susceptible: S
    data_cat[1].push({t: date, y: Math.round(sol_history[i][1]) + Math.round(sol_history[i][2])}); //exposed: E0 + E1
    data_cat[2].push({t: date, y: Math.round(sol_history[i][3])}); //asymptomatic: I0
    data_cat[3].push({t: date, y: Math.round(sol_history[i][5])}); //mild unreported: I1u
    data_cat[4].push({t: date, y: Math.round(sol_history[i][10])}); //recovered unreported: Ru
    data_cat[5].push({t: date, y: Math.round(sol_history[i][9])}); //recovered diagnosed: Rd
    data_cat[6].push({t: date, y: Math.round(sol_history[i][11])}); //fatal: D
    data_cat[7].push({t: date, y: Math.round(sol_history[i][8])}); //critical: I3
    data_cat[8].push({t: date, y: Math.round(sol_history[i][7])}); //severe: I2
    data_cat[9].push({t: date, y: Math.round(sol_history[i][4]) + Math.round(sol_history[i][6])}); //mild diagnosed: I1d + Iq

    let num_confirmed_cases = 0;
    for (let j = 0; j < report_sum_indices.length; ++j)
      num_confirmed_cases += Math.round(sol_history[i][report_sum_indices[j]]);
    data_agg.push({t: date, y: num_confirmed_cases});
  }

  return {total: data_agg, categorized: data_cat};
}

function predictModel(params)
{
  //Rates [1/day]
  const a0 = params.a0;
  const a10 = params.a10;
  const a11 = params.a11;
  const a1 = a10 + a11;
  const g0 = params.g0;
  const g1 = params.g1;
  const p1 = params.p1;
  const g2 = params.g2;
  const p2 = params.p2;
  const g3 = params.g3;
  const mu = params.mu;

  //Initial solution
  let N = params.population;
  let E0_0 = params.E0_0;
  let E1_0 = 0;
  let I0_0 = 0;
  let I1d_0 = 0;
  let I1u_0 = 0;
  let I1q_0 = 0;
  let I2_0 = 0;
  let I3_0 = 0;
  let Rd_0 = params.Rd_0;
  let Ru_0 = 0;
  let D_0 = 0;
  let S_0 = N - E0_0 - E1_0 - I0_0 - I1d_0 - I1u_0 - I1q_0 - I2_0 - I3_0 - Rd_0 - Ru_0 - D_0;

  //Solution vector: [S, E0, E1, I0, I1d, I1u, I1q, I2, I3, Rd, Ru, D]
  let solution_hist = [[S_0, E0_0, E1_0, I0_0, I1d_0, I1u_0, I1q_0, I2_0, I3_0, Rd_0, Ru_0, D_0]];

  const nt = params.T_hist + params.T_pred - 1;
  const nt_sub = 1.0/params.dt;

  for (let i = 0; i < nt; i++)
  {
    let u = solution_hist[i].slice(); //copy last solution vector

    let b1 = params.b1N[i] / N;
    let b2 = params.b2N[i] / N;
    let b3 = params.b3N[i] / N;
    let q_input = params.quarantine_input[i];
    let c = params.diag_frac[i];

    for (let j = 0; j < nt_sub; j++)
    {
      let dS = -(b1*(u[2] + u[3] + u[4] + u[5]) + b2*u[7] + b3*u[8])*u[0];
      let du = [ dS,                                                                                                                  //S
                -dS - a0*u[1],                                                                                                        //E0
                      a0*u[1] - a1*u[2],                                                                                              //E1
                               a10*u[2] - g0*u[3],                                                                                    //I0
                             c*a11*u[2]          - (g1 + p1)*u[4],                                                                    //I1d
                         (1-c)*a11*u[2]                          - (g1 + p1)*u[5],                                                    //I1u
                q_input                                                          - (g1 + p1)*u[6],                                    //I1q
                                                         p1*(u[4]          + u[5]          + u[6]) - (g2 + p2)*u[7],                  //I2
                                                                                                            p2*u[7] - (g3 + mu)*u[8], //I3
                                                         g1*(u[4]                          + u[6])        + g2*u[7]        + g3*u[8], //Rd
                                          g0*u[3] +                       g1*u[5],                                                    //Ru
                                                                                                                             mu*u[8]  //D
               ];

      for (let k = 0; k < u.length; k++)
        u[k] += du[k]*params.dt;
    } //sub-timestepping [hrs]

    if (u[1] < 0.5)
      u[1] = 0.0;

    solution_hist.push(u); //save solution daily
  }

  return solution_hist;
}

function getParameterVector(params)
{
  let n = 2*(params.T_hist - 1) + 6;
  let param_vec = new Array(n).fill(0);
  updateParameterVectorFromStruct(params, param_vec);
  return param_vec;
}

function updateParameterVectorFromStruct(params, param_vec)
{
  let n_b1 = (params.T_hist - 1); //no. of beta_1 values to optimize
  let n_c = (params.T_hist - 1); //no. of c values to optimize

  for (let i = 0; i < n_b1; ++i)
    param_vec[i] = params.b1N[i];

  for (let i = 0; i < n_c; ++i)
    param_vec[n_b1 + i] = params.diag_frac[i];

  param_vec[n_b1 + n_c] = params.E0_0;
  param_vec[n_b1 + n_c + 1] = params.g0;
  param_vec[n_b1 + n_c + 2] = params.g1;
  param_vec[n_b1 + n_c + 3] = params.a0;
  param_vec[n_b1 + n_c + 4] = params.a10;
  param_vec[n_b1 + n_c + 5] = params.a11;
}

function updateParameterStructFromVector(params, param_vec, extrapolate_to_end = false)
{
  let n_b1 = (params.T_hist - 1); //no. of beta_1 values to optimize
  let n_c = (params.T_hist - 1); //no. of c values to optimize

  for (let i = 0; i < n_b1; ++i)
    params.b1N[i] = Math.min(Math.max(param_vec[i], 0.0), 1.0); //beta_1 >= 0

  for (let i = 0; i < n_c; ++i)
    params.diag_frac[i] = Math.min(Math.max(param_vec[n_b1 + i], 0.0), 1.0); // 0 <= c <= 1

  params.E0_0 = Math.max(param_vec[n_b1 + n_c], 1.0);
  params.g0 = Math.min(Math.max(param_vec[n_b1 + n_c + 1], 0.0), 1.0);
  params.g1 = Math.min(Math.max(param_vec[n_b1 + n_c + 2], 0.0), 1.0);
  params.a0 = Math.min(Math.max(param_vec[n_b1 + n_c + 3], 0.0), 1.0);
  params.a10 = Math.min(Math.max(param_vec[n_b1 + n_c + 4], 0.0), 1.0);
  params.a11 = Math.min(Math.max(param_vec[n_b1 + n_c + 5], 0.0), 1.0);

  if (extrapolate_to_end)
  {
    //Extend last beta and c value to future
    for (let i = n_b1; i < params.b1N.length; ++i)
      params.b1N[i] = params.b1N[n_b1-1];

    for (let i = n_c; i < params.diag_frac.length; ++i)
      params.diag_frac[i] = params.diag_frac[n_c-1];
  }
}

function optimizeParameters()
{
  let T_pred_orig = sim_params.T_pred;
  let params = sim_params;
  params.T_pred = 0;
  console.log("Original gamma: " + params.g0 + ", " + params.g1);

  //let params = initializeSimulationParameters(data_real.total.length, 0); //no prediction
  let param_vec = getParameterVector(params);

  let res = getFitResidual(params, param_vec);
  let resnorm = getL2Norm(res);
  let resnorm_init = resnorm;

  const m = res.length;
  const n = param_vec.length;

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
      //Update param_vec based on (scaled) update vector
      for (let i = 0; i < n; ++i)
        param_vec[i] = param_vec0[i] - eta*dparam_vec[i];

      //This checks for the validity of parameters, and may modify values
      updateParameterStructFromVector(params, param_vec);
      updateParameterVectorFromStruct(params, param_vec); //update param_vec in case values were modified

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

  console.log("Final resnorm: " + resnorm + ", rel: " + resnorm/resnorm_init);

  //Copy final solution to global simulation parameters
  updateParameterStructFromVector(sim_params, param_vec, true);
  sim_params.T_pred = T_pred_orig;

  console.log("Final gamma: " + params.g0 + ", " + params.g1);

  updateParameters(true);
  return params;
}

function getFitResidual(params)
{
  const num_eq = 3;
  let sol_hist = predictModel(params);
  let residual = new Array(num_eq*(sol_hist.length-1)).fill(0);

  for (let i = 1; i < sol_hist.length; ++i)
  {
    //Error in number of active patients
    let num_active_pred = sol_hist[i][4] + sol_hist[i][6] + sol_hist[i][7] + sol_hist[i][8]; //I1d + I1q + I2 + I3
    let num_active_true = data_real.categorized[i].y[0] - data_real.categorized[i].y[1] - data_real.categorized[i].y[2];
    residual[num_eq*(i-1)] = num_active_pred - num_active_true;
    // residual[num_eq*(i-1)] = num_active_pred - num_active_true + sol_hist[i][9] - data_real.categorized[i].y[1];
    // residual[num_eq*(i-1)] = num_active_pred + sol_hist[i][9] + sol_hist[i][11] - data_real.categorized[i].y[0];

    //Error in no. of recovered-diagnosed patients
    residual[num_eq*(i-1) + 1] = sol_hist[i][9] - data_real.categorized[i].y[1];

    //Error in no. of fatalities
    residual[num_eq*(i-1) + 2] = sol_hist[i][11] - data_real.categorized[i].y[2];

    if (isNaN(residual[num_eq*(i-1)]) || isNaN(residual[num_eq*(i-1) + 1]) || isNaN(residual[num_eq*(i-1) + 2]))
      console.log("Found NaN");
  }
  return residual;
}

function getFitJacobian(params)
{
  const m = 3*(params.T_hist - 1);

  let param_vec = getParameterVector(params);
  const n = param_vec.length;

  let jac = new Array(m*n).fill(0);

  const delta = 1e-5;

  for (let j = 0; j < n; ++j)
  {
    //Compute finite difference
    param_vec[j] += delta;
    updateParameterStructFromVector(params, param_vec);
    let Rp = getFitResidual(params);

    param_vec[j] -= 2*delta;
    updateParameterStructFromVector(params, param_vec);
    let Rm = getFitResidual(params);

    param_vec[j] += delta;

    for (let i = 0; i < m; ++i)
    {
      jac[i*n + j] = (Rp[i] - Rm[i])/(2*delta);
      if (isNaN(jac[i*n+j]))
        console.log("found nan");
    }
  }

  return jac;
}

function getCurrentPredictionError()
{
  let res_sq = 0.0, res0_sq = 0.0;
  for (let i = 1; i < data_real.total.length; ++i)
  {
    let active_true = data_real.categorized[i].y[0] - data_real.categorized[i].y[1] - data_real.categorized[i].y[2];
    let active_pred = data_predicted.categorized[7][i].y + data_predicted.categorized[8][i].y + data_predicted.categorized[9][i].y;

    let err_a = active_pred - active_true;
    let err_r = data_predicted.categorized[5][i].y - data_real.categorized[i].y[1]; //recovered
    let err_d = data_predicted.categorized[6][i].y - data_real.categorized[i].y[2]; //fatal

    res_sq += err_a*err_a + err_r*err_r + err_d*err_d;
    res0_sq += active_true*active_true + data_real.categorized[i].y[1]*data_real.categorized[i].y[1] + data_real.categorized[i].y[2]*data_real.categorized[i].y[2];
  }
  return Math.sqrt(res_sq)/Math.sqrt(res0_sq);
}

function getL2Norm(vec)
{
  let norm = 0.0;
  for (let i = 0; i < vec.length; ++i)
    norm += vec[i]*vec[i];
  return Math.sqrt(norm);
}
