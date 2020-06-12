// COVID-19 Simulation Tool in JavaScript
// Copyright 2020 Savithru Jayasinghe and Dushan Wadduwage
// Licensed under the MIT License (LICENSE.txt)

'use strict';

// var date_format = 'DD-MM-YYYY';
var date_format = 'YYYY-MM-DD';

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
                 {t: "2020-03-21", y: 6},
                 {t: "2020-05-24", y: 49}, //should be distributed over many days
                 {t: "2020-05-25", y: 41},
                 {t: "2020-05-26", y: 127},
                 {t: "2020-05-27", y: 97},
                 {t: "2020-05-28", y: 35},
                 {t: "2020-05-29", y: 11},
                 {t: "2020-05-30", y: 38},
                 {t: "2020-05-31", y: 9},
                 {t: "2020-06-01", y: 8},
                 {t: "2020-06-02", y: 32},
                 {t: "2020-06-03", y: 35},
                 {t: "2020-06-04", y: 6},
                 {t: "2020-06-05", y: 2},
                 {t: "2020-06-06", y: 4},
                 {t: "2020-06-07", y: 21},
                 {t: "2020-06-08", y: 12},
                 {t: "2020-06-09", y: 1},
                 {t: "2020-06-10", y: 2},]
}

var data_start_dates = {
    "Sri Lanka": "2020-03-01"
};

// var custom_country_data = {
//   "Sri Lanka" : {
//       t_start: [0, 13, 34, 41], //indices to start dates of any interventions
//       b1N: [0.85, 0.14, 0.2, 0.3], //values of b1N for each intervention segment defined in t_start
//       b2N: new Array(4).fill(0), //values of b2N
//       b3N: new Array(4).fill(0),
//       ce: new Array(4).fill(0),
//       c0: new Array(4).fill(0),
//       c1: new Array(4).fill(0.1),
//       c2: new Array(4).fill(1.0),
//       c3: new Array(4).fill(1.0),
//       E0_0: 5, //no. of individuals exposed at start
//       Rd_0: 1, //no. of recovered-diagnosed individuals at start
//   }
// }

// var custom_country_data = {
//   "Sri Lanka" : {
//       t_start: [0, 13, 17, 25, 37, 39, 41, 45, 47, 53, 55], //indices to start dates of any interventions
//       b1N: [0.787, 0.379, 0.093, 0.045, 0.045, 0.045, 0.470, 0.730, 0.730, 0.472, 0.472], //values of b1N for each intervention segment defined in t_start
//       b2N: new Array(11).fill(0), //values of b2N
//       b3N: new Array(11).fill(0),
//       ce: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.114, 0.508],
//       c0: [0, 0, 0, 0, 0, 0, 0, 0.658, 0.801, 0.801, 0.801],
//       c1: [0.029, 0.029, 0.029, 0.029, 0.465, 0.884, 0.884, 0.973, 0.973, 0.973, 0.973],
//       c2: new Array(11).fill(1.0),
//       c3: new Array(11).fill(1.0),
//       E0_0: 5, //no. of individuals exposed at start
//       Rd_0: 1, //no. of recovered-diagnosed individuals at start
//   }
// }

var custom_country_data = {
  "Sri Lanka" : {
              //Mar 1, Mar 15, Mar 25, Apr 10, Apr 15, Apr 25, Apr 30
      t_start: [0, 14, 24, 40, 45, 55, 60], //indices to start dates of any interventions
      b1N: [0.8, 0.1, 0.05, 0.470, 0.6, 0.48, 0.48], //values of b1N for each intervention segment defined in t_start
      b2N: new Array(7).fill(0), //values of b2N
      b3N: new Array(7).fill(0),
      ce: [0, 0, 0, 0, 0, 0.5, 0.5],
      c0: [0, 0, 0, 0, 0.7, 0.7, 0.7],
      c1: [0.03, 0.03, 0.03, 0.965, 0.965, 0.965, 0.1],
      c2: new Array(7).fill(1.0),
      c3: new Array(7).fill(1.0),
      E0_0: 5, //no. of individuals exposed at start
      Rd_0: 1, //no. of recovered-diagnosed individuals at start
  }
}

// var b1_f = 0.3;
// var b1_q90 = 0.14 + 0.1*(b1_f - 0.14);
// var b1_q50 = 0.14 + 0.5*(b1_f - 0.14);
// var b1_cyc = b1_q50;

// var custom_country_data = {
//   "Sri Lanka" : {
//       t_start: [0, 13, 50], //indices to start dates of any interventions
//       b1N: [0.85, 0.14, b1_f], //values of b1N for each intervention segment defined in t_start
//       b2N: new Array(3).fill(0), //values of b2N
//       b3N: new Array(3).fill(0),
//       c1: new Array(3).fill(0.1),
//       E0_0: 5, //no. of individuals exposed at start
//       Rd_0: 1, //no. of recovered-diagnosed individuals at start
//   }
// }

// var custom_country_data = {
//   "Sri Lanka" : {
//       t_start: [0, 13, 50, 84], //indices to start dates of any interventions
//       b1N: [0.85, 0.14, 0.3, 0.14], //values of b1N for each intervention segment defined in t_start
//       b2N: new Array(4).fill(0), //values of b2N
//       b3N: new Array(4).fill(0),
//       c1: new Array(4).fill(0.1),
//       E0_0: 5, //no. of individuals exposed at start
//       Rd_0: 1, //no. of recovered-diagnosed individuals at start
//   }
// }

// var custom_country_data = {
//   "Sri Lanka" : {
//       t_start: [0, 13, 61, 92, 122, 153, 184, 214, 245, 275, 306, 336], //indices to start dates of any interventions
//       b1N: [0.85, 0.14, b1_cyc, 0.14, b1_cyc, 0.14, b1_cyc, 0.14, b1_cyc, 0.14, b1_cyc, 0.14], //values of b1N for each intervention segment defined in t_start
//       b2N: new Array(12).fill(0), //values of b2N
//       b3N: new Array(12).fill(0),
//       c1: new Array(12).fill(0.1),
//       E0_0: 5, //no. of individuals exposed at start
//       Rd_0: 1, //no. of recovered-diagnosed individuals at start
//   }
// }

// var custom_country_data = {
//   "Sri Lanka" : {
//       t_start: [0, 13, 61, 75, 92, 106, 122, 136, 153, 167, 184, 198, 214, 228, 245, 259, 275, 289, 306], //indices to start dates of any interventions
//       b1N: [0.85, 0.14, b1_cyc, 0.14, b1_cyc, 0.14, b1_cyc, 0.14, b1_cyc, 0.14, b1_cyc, 0.14, b1_cyc, 0.14, b1_cyc, 0.14, b1_cyc, 0.14, b1_cyc], //values of b1N for each intervention segment defined in t_start
//       b2N: new Array(19).fill(0), //values of b2N
//       b3N: new Array(19).fill(0),
//       c1: new Array(19).fill(0.1),
//       E0_0: 5, //no. of individuals exposed at start
//       Rd_0: 1, //no. of recovered-diagnosed individuals at start
//   }
// }


//The control parameters will be set to these default values when the user first loads the page.
var default_controls = {
  T_pred: 14,    //prediction length
  b1N: 0.0,     //beta_1 value
  b2N: 0.0,     //beta_2 value
  b3N: 0.0,     //beta_3 value
  ce: 0.0,      //fraction of exposed patients diagnosed daily
  c0: 0.0,      //fraction of asymptomatic patients diagnosed daily
  c1: 0.5,      //fraction of mild patients diagnosed daily
  c2: 1.0,      //fraction of severe patients diagnosed daily
  c3: 1.0,      //fraction of critical patients diagnosed daily
  CFR: 0.02,    //case fatality ratio
  param_error: 0.05 //error in parameters
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

var logarithmic_ticks = {
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

var linear_ticks = {
  min: 0,
  //max: 100000,
  callback: function (value, index, values) {
    if (value == 0) return "0";
    if (value < 1e3) return value.toFixed(0);
    if (value < 1e6) return (value/1e3).toFixed(1) + "k";
    if (value < 1e9) return (value/1e6).toFixed(1) + "M";
    if (value < 1e12) return (value/1e9).toFixed(1) + "B";
    return null;
  }
};


window.onload = function()
{
  //TODO: Temporary code for generating spikes in beta-1
  // for (let t = 61; t < 350; t += 14)
  // {
  //   custom_country_data["Sri Lanka"].t_start.push(t);
  //   custom_country_data["Sri Lanka"].t_start.push(t+1);
  //   custom_country_data["Sri Lanka"].b1N.push(0.60);
  //   custom_country_data["Sri Lanka"].b1N.push(0.14);
  // }
  // custom_country_data["Sri Lanka"].b2N = new Array(custom_country_data["Sri Lanka"].t_start.length).fill(0);
  // custom_country_data["Sri Lanka"].b3N = new Array(custom_country_data["Sri Lanka"].t_start.length).fill(0);
  // custom_country_data["Sri Lanka"].c1 = new Array(custom_country_data["Sri Lanka"].t_start.length).fill(0.1);

  generateCountryDropDown();
  changeCountry(active_country);

  //Update UI controls to match default values
  document.getElementById("slider_finalT").value = default_controls.T_pred;
  document.getElementById("slider_finalT_value").innerHTML = default_controls.T_pred;

  document.getElementById("slider_param_CFR").value = default_controls.CFR * 100;
  document.getElementById("slider_param_CFR_value").innerHTML = (default_controls.CFR * 100).toFixed(1);

  document.getElementById("slider_param_error").value = default_controls.param_error * 100;
  document.getElementById("slider_param_error_value").innerHTML = (default_controls.param_error * 100).toFixed(1);

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
  document.getElementById("R0_value").innerHTML = getR0(sim_params, 0).toFixed(1);

  if (main_chart)
    refreshAllChartData();
  else {
    setupChart();
    setupControlChart();
  }
}

function getCountryData(country_name)
{
  let data_array = world_data[country_name];

  let start_date = data_start_dates[country_name];
  if (start_date)
    start_date = moment(start_date, date_format);

  let data_cat = [], data_total = [], data_fatal = [];
  if (start_date)
  {
    for (let data of data_array)
    {
      let data_t = moment(data.date, date_format);
      if (start_date <= data_t)
      {
        data_total.push({t: data_t, y: data.confirmed});
        data_cat.push({t: data_t, y: [data.confirmed, data.recovered, data.deaths, 0]});
        data_fatal.push({t: data_t, y: data.deaths});
      }
    }
  }
  else
  {
    let first_time = true;
    for (let data of data_array)
    {
      let data_t = moment(data.date, date_format);
      if (data.confirmed + data.recovered + data.deaths > 0)
      {
        if (first_time) //Add buffer days before first confirmed case
        {
          for (let i = 5; i > 0; --i)
          {
            let t_tmp = data_t.clone().subtract(i, 'days');
            data_total.push({t: t_tmp, y: 0});
            data_cat.push({t: t_tmp, y: [0, 0, 0, 0]});
            data_fatal.push({t: t_tmp, y: 0});
          }
          first_time = false;
        }
        data_total.push({t: data_t, y: data.confirmed});
        data_cat.push({t: data_t, y: [data.confirmed, data.recovered, data.deaths, 0]});
        data_fatal.push({t: data_t, y: data.deaths});
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

  return {total: data_total, categorized: data_cat, fatal: data_fatal};
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
        params.ce[j] = data.ce[i];
        params.c0[j] = data.c0[i];
        params.c1[j] = data.c1[i];
        params.c2[j] = data.c2[i];
        params.c3[j] = data.c3[i];
      }
    }

    //Update initial E0 and Rd
    if (data.E0_0)
      params.E0_0 = data.E0_0;

    if (data.Rd_0)
      params.Rd_0 = data.Rd_0;
  }

  //Update quarantine input data
  for (let i = 0; i < data_real.categorized.length; ++i)
    params.quarantine_input[i] = data_real.categorized[i].y[3];

  //Update population
  let N = world_population[country_name];
  if (N)
    params.population = N;
  else
    console.log("Population data not found for " + country_name + "!");
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
  else {
    yaxis_config.type = 'linear';
    yaxis_config.ticks = linear_ticks;
  }

  const bar_width_frac = 1.0;
  const cat_width_frac = 0.9;

	let main_chart_config = {
			data: {
				datasets: [
				// {
				// 	label: 'Susceptible',
				// 	backgroundColor: '#eeeeec',
				// 	data: data_predicted.categorized[0],
				// 	type: 'bar',
				// 	stack: 'stack0',
				// 	hidden: true,
				// 	barPercentage: bar_width_frac,
				// 	categoryPercentage: cat_width_frac,
				// 	order: 13
				// },
				{
					label: 'Exposed',
					backgroundColor: 'rgba(255, 220, 160, 0.75)',
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 10
				},
        {
          label: 'Asymptomatic',
          backgroundColor: 'rgba(180, 240, 255, 0.75)',
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 9
        },
        {
          label: 'Recovered',
          backgroundColor: 'rgba(130, 210, 50, 0.75)', //#73d216',
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 8
        },
        {
          label: 'Fatal',
          backgroundColor: 'rgba(10, 10, 10, 0.75)',
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 7
        },
        {
          label: 'Mildly infected',
          backgroundColor: 'rgba(255, 200, 0, 0.75)',
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 6
        },
        {
          label: 'Severely infected',
          backgroundColor: 'rgba(240, 150, 40, 0.75)',
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 5
        },
        {
          label: 'Critically infected',
          backgroundColor: 'rgba(200, 0, 0, 0.75)',
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 4
        },
				{
					label: 'Actual reported',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#3465a4',
					type: 'line',
					fill: true,
					borderWidth: 2,
          pointRadius: 2,
					order: 0
				},
				{
					label: 'Predicted reported',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#729fcf',
					borderDash: [5, 5],
					type: 'line',
					fill: true,
					borderWidth: 2,
          pointRadius: 2,
					order: 1
				},
        {
          label: 'Actual deaths',
          backgroundColor: 'rgba(1,1,1,0)',
          borderColor: 'rgb(10, 10, 10)',
          hidden: true,
          type: 'line',
          fill: true,
          borderWidth: 2,
          pointRadius: 2,
          order: 2
        },
        {
          label: 'Predicted deaths',
          backgroundColor: 'rgba(1,1,1,0)',
          borderColor: 'rgb(80, 80, 80)',
          borderDash: [5, 5],
          hidden: true,
          type: 'line',
          fill: true,
          borderWidth: 2,
          pointRadius: 2,
          order: 3
        },
        {
          label: 'Predicted lower',
          backgroundColor: 'rgba(200,200,200,0.4)',
          borderColor: 'rgba(100,100,100, 0.4)',
          type: 'line',
          fill: '+1',
          borderWidth: 1,
          pointRadius: 0,
          order: 11
        },
        {
          label: 'Predicted upper',
          backgroundColor: 'rgba(200,200,200,0.4)',
          borderColor: 'rgba(100,100,100, 0.4)',
          type: 'line',
          fill: '-1',
          borderWidth: 1,
          pointRadius: 0,
          order: 12
        },
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
                    updateLegend(tooltipItem.index);
                    document.getElementById("R0_value").innerHTML = getR0(sim_params, tooltipItem.index).toFixed(1);

                    if (tooltipItem.datasetIndex >= 7)
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

  refreshMainChartData();
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
        labelString: 'Beta1',
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
          label: 'Beta1',
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
            if (el[0])
            {
              control_chart_canvas.style.cursor = "pointer";
              document.getElementById("R0_value").innerHTML = getR0(sim_params, el[0]._index).toFixed(1);
            }
            else {
              control_chart_canvas.style.cursor = "default";
            }
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
          if (Math.abs(datavec[i] - yval_old) > 2e-3)
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
    let array_name = document.querySelector('input[name="plot_data_type"]:checked').value; //returns cat_diag or cat_sum

    let n_cat = data_predicted[array_name].length;
    for (let i = 0; i < n_cat; ++i)
      main_chart.data.datasets[i].data = data_predicted[array_name][i];

    main_chart.data.datasets[n_cat].data = data_real.total;
    main_chart.data.datasets[n_cat+1].data = data_predicted.total;

    main_chart.data.datasets[n_cat+2].data = data_real.fatal;
    main_chart.data.datasets[n_cat+3].data = data_predicted[array_name][3];

    main_chart.data.datasets[n_cat+4].data = data_predicted.lower;
    main_chart.data.datasets[n_cat+5].data = data_predicted.upper;

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
    "ce": "c_e",
    "c0": "c_0",
    "c1": "c_1",
    "c2": "c_2",
    "c3": "c_3"
  };

  const param_to_titlestring = {
    "b1N": "Transmission rate of unreported exposed, asymptomatic, or mildly-infected individuals",
    "b2N": "Transmission rate of unreported severely-infected individuals",
    "b3N": "Transmission rate of unreported critically-infected individuals",
    "ce": "Fraction of exposed cases diagnosed daily",
    "c0": "Fraction of asymptomatic cases diagnosed daily",
    "c1": "Fraction of mildly-infected cases diagnosed daily",
    "c2": "Fraction of severely-infected cases diagnosed daily",
    "c3": "Fraction of critically-infected cases diagnosed daily"
  };

  if (control_chart)
  {
    control_chart.data.datasets[0].data = getControlChartData();
    control_chart.data.datasets[0].label = param_to_labelstring[active_control_parameter];
    control_chart.options.scales.yAxes[0].scaleLabel.labelString = param_to_labelstring[active_control_parameter];
    control_chart.update();
    delete control_chart.$zoom._originalOptions[control_chart.options.scales.xAxes[0].id].time.min;
    delete control_chart.$zoom._originalOptions[control_chart.options.scales.xAxes[0].id].time.max;
    document.getElementById("control_chart_title").innerHTML = param_to_titlestring[active_control_parameter];
  }
}

function setLogYAxis(is_log)
{
  if (is_log)
  {
    main_chart.options.scales.yAxes[0].type = 'logarithmic';
    main_chart.options.scales.yAxes[0].ticks = logarithmic_ticks;
  }
  else
  {
    main_chart.options.scales.yAxes[0].type = 'linear';
    main_chart.options.scales.yAxes[0].ticks = linear_ticks;
  }
  main_chart.update();
}

function toggleDatasets()
{
  let calibration_mode = document.getElementById("check_calibration_mode").checked;

  for (let i = 0; i < 7; ++i)
    main_chart.data.datasets[i].hidden = calibration_mode;

  for (let i = 9; i < 11; ++i)
    main_chart.data.datasets[i].hidden = !calibration_mode; //actual deaths, predicted deaths

  for (let i = 11; i < 13; ++i)
    main_chart.data.datasets[i].hidden = calibration_mode;

  main_chart.update();
}

function resetZoom()
{
  main_chart.resetZoom();
  control_chart.resetZoom();
}

function resetParameters()
{
  sim_params = initializeSimulationParameters(data_real.total.length, default_controls.T_pred);
  updateParameters(true);
}

function updateLegend(day = last_active_tooltip_day)
{
  if (day < 0 || day >= data_predicted.total.length)
    return;

  last_active_tooltip_day = day;

  document.getElementById("legend_date").innerHTML = data_predicted.total[day].t.format("MMM-DD-YYYY");

  let num_total = 0;
  for (let i = 0; i < data_predicted.cat_diag.length; ++i)
  {
    document.getElementById("legend_predu" + i).innerHTML = formatNumber(data_predicted.cat_sum[i][day].y - data_predicted.cat_diag[i][day].y);
    document.getElementById("legend_pred" + i).innerHTML = formatNumber(data_predicted.cat_diag[i][day].y);
    num_total += data_predicted.cat_sum[i][day].y;
  }

  let num_infected_diag = data_predicted.cat_diag[4][day].y + data_predicted.cat_diag[5][day].y + data_predicted.cat_diag[6][day].y;
  let num_infected_sum = data_predicted.cat_sum[4][day].y + data_predicted.cat_sum[5][day].y + data_predicted.cat_sum[6][day].y;
  document.getElementById("legend_predu_infected").innerHTML = formatNumber(num_infected_sum - num_infected_diag);
  document.getElementById("legend_pred_infected").innerHTML = formatNumber(num_infected_diag);

  document.getElementById("legend_predu_total").innerHTML = formatNumber(num_total - data_predicted.total[day].y);
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
  let total_length = hist_length + 280;

  let params = {
    T_hist: hist_length,
    T_pred: pred_length,
    dt: 1.0/24.0,                                                       //timestep size [days]
    param_error: default_controls.param_error,                          //assumed error in rate parameters
    b1N: new Array(total_length).fill(default_controls.b1N),            //transmission rate from mild to susceptible
    b2N: new Array(total_length).fill(default_controls.b2N),            //transmission rate from severe to susceptible
    b3N: new Array(total_length).fill(default_controls.b3N),            //transmission rate from critical to susceptible
    quarantine_input: new Array(total_length).fill(0.0),                //no. of patients added directly to quarantine
    ce: new Array(total_length).fill(default_controls.ce),              //fraction of exposed patients diagnosed daily
    c0: new Array(total_length).fill(default_controls.c0),              //fraction of asymptomatic patients diagnosed daily
    c1: new Array(total_length).fill(default_controls.c1),              //fraction of mild patients diagnosed daily
    c2: new Array(total_length).fill(default_controls.c2),              //fraction of severe patients diagnosed daily
    c3: new Array(total_length).fill(default_controls.c3),              //fraction of critical patients diagnosed daily
    T_incub0: 3,                                                        //periods [days]
    T_incub1: 2,
    T_asympt: 6,
    T_mild  : 6,
    T_severe: 4,
    T_icu   : 10,
    f: 0.3,                                                             //exposed to asymptomatic probability
    frac_recover_I1: 0.80,                                              //fraction of cases that recover from mild-infected stage I1
    frac_recover_I2: 0.75,                                              //fraction of cases that recover from severe-infected stage I2
    frac_recover_I3: 0.60,                                              //fraction of cases that recover from critical-infected stage I3
    population: 1E7,                                                    //population of country
    E0_0: 5,                                                            //number of non-infectious exposed individuals at start
    Rd_0: 0,                                                            //number of recovered-diagnosed individuals at start
  }

  setRateParameters(params);

  //Modify any parameters that are specific to the currently active country.
  customizeParametersByCountry(active_country, params);

  return params;
}

function setRateParameters(params)
{
  //Probabilities
  const prob_E0_E1 = 1;  //non-infectious exposed to infectious exposed
  const prob_E1_I0 = params.f;  //exposed to asymptomatic
  const prob_E1_I1 = 1 - prob_E1_I0;  //exposed to mild
  const prob_I0_R  = 1;
  const prob_I1_R  = params.frac_recover_I1;  //mild to recovered //0.8
  const prob_I1_I2 = 1 - prob_I1_R; //mild to severe  //0.2
  const prob_I2_R  = params.frac_recover_I2; //severe to recovered  //0.75
  const prob_I2_I3 = 1 - prob_I2_R; //severe to critical //0.25
  const prob_I3_R  = params.frac_recover_I3;; //critical to recovered //0.6
  const prob_I3_D  = 1 - prob_I3_R; //critical to dead //0.4


  //set rate parameters [1/day]
  params.a0  = (1/params.T_incub0) * prob_E0_E1;
  params.a10 = (1/params.T_incub1) * prob_E1_I0;
  params.a11 = (1/params.T_incub1) * prob_E1_I1;
  params.g0  = (1/params.T_asympt) * prob_I0_R;
  params.g1  = (1/params.T_mild)   * prob_I1_R;
  params.p1  = (1/params.T_mild)   * prob_I1_I2;
  params.g2  = (1/params.T_severe) * prob_I2_R;
  params.p2  = (1/params.T_severe) * prob_I2_I3;
  params.g3  = (1/params.T_icu)    * prob_I3_R;
  params.mu  = (1/params.T_icu)    * prob_I3_D;
}

//Computes the basic reproduction number (R0) from the model parameters
function getR0(params, ind)
{
  let b1N = params.b1N[ind];
  let b2N = params.b2N[ind];
  let b3N = params.b3N[ind];
  let a1 = params.a10 + params.a11;
  let tmp = b1N + params.p1*(b2N + b3N*params.p2/(params.mu + params.g3)) / (params.p2 + params.g2);
  return (b1N/a1 + b1N*params.f/params.g0 + (1-params.f)/(params.p1 + params.g1)*tmp);
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

  let slider_param_CFR = document.getElementById("slider_param_CFR");
  if (slider_param_CFR)
  {
    let CFR = Number(slider_param_CFR.value) / 100.0;
    let frac_recover_I3 = 1.0 - CFR / ((1.0 - sim_params.frac_recover_I1)*(1.0 - sim_params.frac_recover_I2));
    if (sim_params.frac_recover_I3 != frac_recover_I3)
    {
      sim_params.frac_recover_I3 = frac_recover_I3;
      setRateParameters(sim_params);
      requires_update = true;
      document.getElementById("slider_param_CFR_value").innerHTML = (CFR*100).toFixed(1);
    }
  }

  let slider_param_error = document.getElementById("slider_param_error");
  if (slider_param_error)
  {
    let val = Number(slider_param_error.value) / 100.0;
    if (sim_params.param_error != val)
    {
      sim_params.param_error = val;
      requires_update = true;
      document.getElementById("slider_param_error_value").innerHTML = (val*100).toFixed(1);
    }
  }

  if (requires_update)
  {
    data_predicted = getPredictionData(data_real.total[0].t);
    refreshAllChartData();
    document.getElementById("prediction_error").innerHTML = getCurrentPredictionError().toFixed(3);
  }
}

function getErrorData()
{
  //Save off original parameters
  const total_length = sim_params.b1N.length;
  let b1N = sim_params.b1N.slice();
  let b2N = sim_params.b2N.slice();
  let b3N = sim_params.b3N.slice();

  const names = ["a0", "a10", "a11", "g0", "g1", "p1", "g2", "p2", "g3", "mu"];

  let params_orig = {};
  for (let name of names)
    params_orig[name] = sim_params[name];

  const f_upper = 1 + sim_params.param_error;
  const f_lower = 1 - sim_params.param_error;

  for (let i = 0; i < total_length; ++i)
  {
    sim_params.b1N[i] = b1N[i] * f_lower;
    sim_params.b2N[i] = b2N[i] * f_lower;
    sim_params.b3N[i] = b3N[i] * f_lower;
  }
  for (let name of names)
    sim_params[name] = params_orig[name] * f_lower;

  let sol_history_lower = predictModel(sim_params);

  for (let i = 0; i < total_length; ++i)
  {
    sim_params.b1N[i] = b1N[i] * f_upper;
    sim_params.b2N[i] = b2N[i] * f_upper;
    sim_params.b3N[i] = b3N[i] * f_upper;
  }
  for (let name of names)
    sim_params[name] = params_orig[name] * f_upper;

  let sol_history_upper = predictModel(sim_params);

  //Restore original parameters
  for (let i = 0; i < total_length; ++i)
  {
    sim_params.b1N[i] = b1N[i];
    sim_params.b2N[i] = b2N[i];
    sim_params.b3N[i] = b3N[i];
  }
  for (let name of names)
    sim_params[name] = params_orig[name];

  return {lower: sol_history_lower, upper: sol_history_upper};
}

function getPredictionData(start_date)
{
  let pop_hist = predictModel(sim_params); //population history
  let error_data = getErrorData();

  let data_agg = [], data_lower = [], data_upper = [];
  let data_cat_diag = new Array(7);
  let data_cat_sum = new Array(data_cat_diag.length);
  for (let i = 0; i < data_cat_diag.length; ++i)
  {
    data_cat_diag[i] = new Array();
    data_cat_sum[i] = new Array();
  }

  // const report_sum_indices = [4, 6, 7, 8, 9, 11]; //I1d + I1q + I2 + I3 + Rd + D

  for (let i = 0; i < pop_hist.length; i++)
  {
    let date = start_date.clone().add(i,'days');

    //Accumulate data into categories for plotting
    data_cat_diag[0].push({t: date, y: Math.round(pop_hist[i].E1[1])}); //exposed: E1d
    data_cat_diag[1].push({t: date, y: Math.round(pop_hist[i].I0[1])}); //asymptomatic: I0d
    data_cat_diag[2].push({t: date, y: Math.round(pop_hist[i].R[1])}); //recovered: Rd
    data_cat_diag[3].push({t: date, y: Math.round(pop_hist[i].D[1])}); //fatal: D
    data_cat_diag[4].push({t: date, y: Math.round(pop_hist[i].I1[1])}); //mild: I1d
    data_cat_diag[5].push({t: date, y: Math.round(pop_hist[i].I2[1])}); //severe: I2d
    data_cat_diag[6].push({t: date, y: Math.round(pop_hist[i].I3[1])}); //critical: I3d

    data_cat_sum[0].push({t: date, y: pop_hist[i].getNumExposed()}); //exposed: E0 + E1u + E1d
    data_cat_sum[1].push({t: date, y: pop_hist[i].getNumAsymptomatic()}); //asymptomatic: I0u + I0d
    data_cat_sum[2].push({t: date, y: pop_hist[i].getNumRecovered()}); //recovered: Ru + Rd
    data_cat_sum[3].push({t: date, y: pop_hist[i].getNumFatal()}); //fatal: Du + Dd
    data_cat_sum[4].push({t: date, y: pop_hist[i].getNumMild()}); //mild: I1u + I1d
    data_cat_sum[5].push({t: date, y: pop_hist[i].getNumSevere()}); //severe: I2u + I2d
    data_cat_sum[6].push({t: date, y: pop_hist[i].getNumCritical()}); //critical: I3u + I3d

    let num_diag = pop_hist[i].getNumDiagnosed();
    let num_diag_lower = error_data.lower[i].getNumDiagnosed();
    let num_diag_upper = error_data.upper[i].getNumDiagnosed();

    data_agg.push({t: date, y: num_diag});
    data_lower.push({t: date, y: Math.min(num_diag, num_diag_lower, num_diag_upper)});
    data_upper.push({t: date, y: Math.max(num_diag, num_diag_lower, num_diag_upper)});
  }

  return {total: data_agg, cat_diag: data_cat_diag, cat_sum: data_cat_sum, lower: data_lower, upper: data_upper};
}

function predictModel(params)
{
  //Create initial population object
  let pop0 = new Population(params.population);
  pop0.E0 = params.E0_0;
  pop0.R[1] = params.Rd_0; //recovered-diagnosed
  pop0.S = pop0.N - pop0.E0 - pop0.R[1];

  let population_hist = [pop0];

  const nt = params.T_hist + params.T_pred - 1;
  const nt_sub = 1.0/params.dt;

  for (let t = 0; t < nt; t++)
  {
    let pop_new = population_hist[t].clone();

    for (let j = 0; j < nt_sub; j++) //sub-timestepping [hrs]
      pop_new.evolve(params, t);

    pop_new.report(params, t); //report once daily

    population_hist.push(pop_new); //save solution daily
  }

  return population_hist;
}
