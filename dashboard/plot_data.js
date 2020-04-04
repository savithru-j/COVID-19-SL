// COVID-19 Simulation Tool in JavaScript
// Copyright 2020 Savithru Jayasinghe and Dushan Wadduwage
// Licensed under the MIT License (LICENSE.txt)

'use strict';
$(function () {

    var date_format = 'YYYY-MM-DD';
    var data_raw_SL = {};
    var data_SL = {};
    var sim_params = {};
    var data_predicted = {};
    var chart = [];
    var chart_config = [];
    var control_chart = [];

    //The control parameters will be set to these default values when the user first loads the page.
    var default_controls = {
        T_pred: 7,        //prediction length
        b1N_init: 0.92,   //initial beta_1 value
        b2N_init: 0.00,   //initial beta_2 value
        b3N_init: 0.00,   //initial beta_3 value
        b1N_T0: 0.05,     //beta_1 after intervention 0
        b1N_T1: 0.05,     //beta_1 after intervention 1
        T0: 13,           //time of intervention 0
        diag_frac: 0.12
    };
    var last_active_tooltip_day = 0;
    var startDataDate = {
        "Sri Lanka":"2020-02-29",
        "Singapore":"2020-02-29"
    };


    generateCountryDropDown();
    showChartForCountry("Singapore");

    function showChartForCountry(c) {
        chart=[];
        chart_config=[];
        selectDropDown(c);
        data_raw_SL = getRaw(c);
        data_SL=getDataSriLanka();
        sim_params = initializeSimulationParameters(data_SL.length, default_controls.T_pred);
        data_predicted = getPredictionData(data_SL[0].t);

        //Update UI controls to match default values
        document.getElementById("slider_finalT").value = default_controls.T_pred;
        document.getElementById("slider_finalT_value").innerHTML = default_controls.T_pred;

        document.getElementById("slider_b1").value = default_controls.b1N_init;
        document.getElementById("slider_b1_value").innerHTML = default_controls.b1N_init.toFixed(2);

        document.getElementById("slider_b2").value = default_controls.b2N_init;
        document.getElementById("slider_b2_value").innerHTML = default_controls.b2N_init.toFixed(2);

        document.getElementById("slider_b3").value = default_controls.b3N_init;
        document.getElementById("slider_b3_value").innerHTML = default_controls.b3N_init.toFixed(2);

        document.getElementById("slider_c").value = default_controls.diag_frac;
        document.getElementById("slider_c_value").innerHTML = default_controls.diag_frac.toFixed(2);

        let slider_interv0_T = document.getElementById('slider_interv0_T');
        let slider_interv1_T = document.getElementById('slider_interv1_T');
        slider_interv0_T.max = sim_params.T_hist + sim_params.T_pred;
        slider_interv1_T.max = sim_params.T_hist + sim_params.T_pred;
        slider_interv0_T.value = default_controls.T0;
        slider_interv1_T.value = slider_interv0_T.max;
        document.getElementById('slider_interv0_T_value').innerHTML = slider_interv0_T.value;
        document.getElementById('slider_interv1_T_value').innerHTML = slider_interv1_T.value;

        document.getElementById("slider_interv0_b1").value = default_controls.b1N_T0;
        document.getElementById("slider_interv0_b1_value").innerHTML = default_controls.b1N_T0.toFixed(2);

        document.getElementById("slider_interv1_b1").value = default_controls.b1N_T1;
        document.getElementById("slider_interv1_b1_value").innerHTML = default_controls.b1N_T1.toFixed(2);

        updateChart();
        updateLegend();
        updateParameters();
        updateControlChart();
    }

    $("#countries").on('change', function (e) {
        var valueSelected = this.value;
        console.log(valueSelected);
        showChartForCountry(valueSelected);
    });

    function getDataSriLanka()
    {
      let data = [];
      for (let i = 0; i < data_raw_SL.length; ++i)
        data.push({t: data_raw_SL[i].t, y: data_raw_SL[i].y[0]});
      return data;
    }

    function generateCountryDropDown() {
        $.each(countries, function( index, item ) {
            var o = new Option(index, index);
            $(o).html(index);
            $("#countries").append(o);
        });
    }

    function selectDropDown(c) {
        $("#countries").val(c);
    }

    //data location
    //https://github.com/pomber/covid19
    //https://pomber.github.io/covid19/timeseries.json
    function getRaw(c){
        var countryData = countries[c];
        var d=[];
        $.each(countryData, function( index, item ) {
            if(startDataDate[c]){
                if(new Date(startDataDate[c])<new Date(item.date)){
                    d.push({t: moment(item.date, date_format), y: [item.confirmed, item.recovered, item.deaths, 0]})
                }
            }else{
                d.push({t: moment(item.date, date_format), y: [item.confirmed, item.recovered, item.deaths, 0]})
            }
        });
        return d;
    }

    $('#check_log_y').change(function() {
        setLogYAxis(this.checked);
    });
    $("#reset_zoom").click(function () {
        resetZoom();
    });
    $("#slider_finalT,#slider_b1,#slider_b2,#slider_b3,#slider_interv0_b1,#slider_interv1_b1,#slider_c").mousemove(function( event ) {
        updateParameters();
    });
    $("#slider_interv0_T,#slider_interv1_T").mousemove(function( event ) {
        updateParameters(true);
    });
    $('#check_interv0').change(function() {
        updateInterventions(0);
    });
    $('#check_interv1').change(function() {
        updateInterventions(1);
    });

    function formatNumber(num) {
        return num.toString().replace(/(\d)(?=(\d{3})+(?!\d))/g, '$1,')
    }

    function updateChart()
    {
        var canvas = document.createElement('canvas');
        canvas.className="full";
        $("#canvas_wrapper").html(canvas);

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
          label: 'Mildly infected - hidden',
          backgroundColor: 'rgba(255, 200, 0, 0.5)',
          data: data_predicted.categorized[3],
          type: 'bar',
          stack: 'stack0',
          barPercentage: bar_width_frac,
          categoryPercentage: cat_width_frac,
          order: 8
        },
        {
          label: 'Recovered - hidden',
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
					label: 'Mildly infected',
					backgroundColor: 'rgba(255, 200, 0, 0.75)',
					data: data_predicted.categorized[6],
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 5
				},
				{
					label: 'Severely infected',
					backgroundColor: 'rgba(240, 150, 40, 0.75)',
					data: data_predicted.categorized[7],
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 4
				},
				{
					label: 'Critically infected',
					backgroundColor: 'rgba(200, 0, 0, 0.75)',
					data: data_predicted.categorized[8],
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 3
				},
				{
					label: 'Fatal',
					backgroundColor: 'rgba(10, 10, 10, 0.75)',
					data: data_predicted.categorized[9],
					type: 'bar',
					stack: 'stack0',
					barPercentage: bar_width_frac,
					categoryPercentage: cat_width_frac,
					order: 2
				},
				{
					label: 'actual diagnosed',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#3465a4',
					data: data_SL,
					type: 'line',
					fill: true,
					borderWidth: 2,
					order: 0
				},
				{
					label: 'predicted diagnosed',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#729fcf',
					borderDash: [5, 5],
					data: data_predicted.aggregated,
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
            boxWidth: 10,
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
                value: data_predicted.aggregated[default_controls.T0-1].t,
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
                        onPan: function () {
                let minVal = chart.scales['x-axis-0']._table[0].time;
                let maxVal = chart.scales['x-axis-0']._table[1].time;
                control_chart.options.scales.xAxes[0].time.min = minVal;
                control_chart.options.scales.xAxes[0].time.max = maxVal;
                control_chart.update();
              }
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
                        onZoom: function () {
                let minVal = chart.scales['x-axis-0']._table[0].time;
                let maxVal = chart.scales['x-axis-0']._table[1].time;
                control_chart.options.scales.xAxes[0].time.min = minVal;
                control_chart.options.scales.xAxes[0].time.max = maxVal;
                control_chart.update();
              }
                    }
                }
            }

                }
            };

      chart = new Chart(ctx, chart_config);
    };

    function updateControlChart()
    {
      let canvas = document.getElementById('control_chart_canvas');
      //canvas.width = 0.7*window.innerWidth;
      canvas.height = 0.2*window.innerHeight;
      let ctx = canvas.getContext('2d');

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
          }
        };

      let control_chart_config = {
          data: {
            datasets: [{
              label: 'beta_1',
              backgroundColor: 'rgba(1,1,1,0)',
              borderColor: 'rgb(90, 110, 150)',
              data: getBeta1Data(),
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
                  onPan: function () {
                    let minVal = control_chart.scales['x-axis-0']._table[0].time;
                    let maxVal = control_chart.scales['x-axis-0']._table[1].time;
                    chart.options.scales.xAxes[0].time.min = minVal;
                    chart.options.scales.xAxes[0].time.max = maxVal;
                    chart.update();
                  }
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
                  onZoom: function () {
                    let minVal = control_chart.scales['x-axis-0']._table[0].time;
                    let maxVal = control_chart.scales['x-axis-0']._table[1].time;
                    chart.options.scales.xAxes[0].time.min = minVal;
                    chart.options.scales.xAxes[0].time.max = maxVal;
                    chart.update();
                  }
                },
              }
            }

          }
        };
      control_chart = new Chart(ctx, control_chart_config);
    }

    function getBeta1Data()
    {
      let data = [];
      for (let i = 0; i < data_predicted.aggregated.length; ++i)
        data.push({t:  data_predicted.aggregated[i].t, y: sim_params.b1N[i]});
      return data;
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

    function resetZoom()
    {
      chart.resetZoom();
      control_chart.resetZoom();
    }

    function updateLegend(day = last_active_tooltip_day)
    {
      if (day < 0 || day >= data_predicted.aggregated.length)
        return;

      last_active_tooltip_day = day;

      document.getElementById("legend_date").innerHTML = data_predicted.aggregated[day].t.format("MMM-DD-YYYY");

      for (let i = 1; i < data_predicted.categorized.length; ++i)
        document.getElementById("legend_pred" + i).innerHTML = formatNumber(data_predicted.categorized[i][day].y);
      document.getElementById("legend_pred_infected").innerHTML = formatNumber(data_predicted.categorized[6][day].y + data_predicted.categorized[7][day].y + data_predicted.categorized[8][day].y);
      document.getElementById("legend_pred_total").innerHTML = formatNumber(data_predicted.aggregated[day].y);

      let true_data = ['-', '-', '-', '-'];
      if (day < data_raw_SL.length)
      {
        true_data[0] = formatNumber(data_raw_SL[day].y[0] - data_raw_SL[day].y[1] - data_raw_SL[day].y[2]);
        true_data[1] = formatNumber(data_raw_SL[day].y[1]);
        true_data[2] = formatNumber(data_raw_SL[day].y[2]);
        true_data[3] = formatNumber(data_raw_SL[day].y[0]);
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

      let q_input = new Array(total_length).fill(0.0);
      for (let i = 0; i < data_raw_SL.length; ++i)
        q_input[i] = data_raw_SL[i].y[3];

      let params = {
        T_hist: hist_length,
        T_pred: pred_length,
        dt: 0.5/24.0,                                                       //timestep size [days]
        b1N: new Array(total_length).fill(default_controls.b1N_init),       //transmission rate from mild to susceptible
        b2N: new Array(total_length).fill(default_controls.b2N_init),       //transmission rate from severe to susceptible
        b3N: new Array(total_length).fill(default_controls.b3N_init),       //transmission rate from critical to susceptible
        quarantine_input: q_input,                                          //no. of patients added directly to quarantine
        diag_frac: new Array(total_length).fill(default_controls.diag_frac),//fraction of I1 patients that are diagnosed
        E0_0: 5,                                                            //number of non-infectious exposed individuals at start
      }

      for (let i = default_controls.T0; i < total_length; ++i)
        params.b1N[i] = default_controls.b1N_T0; //To simulate curfew

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

      let slider_element_ids = ["slider_b2", "slider_b3", "slider_c"];
      let param_arrays = [sim_params.b2N, sim_params.b3N, sim_params.diag_frac];

      for (let i = 0; i < 3; ++i)
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
          document.getElementById("slider_finalT_value").innerHTML = val;
          slider_interv0_T.max = sim_params.T_hist + sim_params.T_pred;
          slider_interv1_T.max = sim_params.T_hist + sim_params.T_pred;
        }
      }

      if (requires_update)
      {
        data_predicted = getPredictionData(data_SL[0].t);
        for (let i = 0; i < data_predicted.categorized.length; ++i)
          chart_config.data.datasets[i].data = data_predicted.categorized[i];

        chart_config.data.datasets[11].data = data_predicted.aggregated;

        for (let i = 0; i < 2; ++i)
        {
          if (interv_params[i].t != Infinity)
            chart.annotation.elements["vertline_T" + i].options.value = data_predicted.aggregated[interv_params[i].t-1].t;
        }
        chart.update();
        updateLegend();

          control_chart.config.data.datasets[0].data = getBeta1Data();
          control_chart.update();
      }
    }

    function getPredictionData(start_date)
    {
      let sol_history = predictModel(sim_params);

      let data_agg = [];
      let data_cat = new Array(10);
      for (let i = 0; i < data_cat.length; ++i)
        data_cat[i] = new Array();

      const report_sum_indices = [5, 6, 7, 8, 10]; //omits c*I1, which will be added inside the loop

      for (let i = 0; i < sol_history.length; i++)
      {
        let date = start_date.clone().add(i,'days');
        let c = sim_params.diag_frac[i];

        //Accumulate data into categories for plotting
        data_cat[0].push({t: date, y: Math.round(sol_history[i][0])}); //S
        data_cat[1].push({t: date, y: Math.round(sol_history[i][1]) + Math.round(sol_history[i][2])}); //E0 + E1
        data_cat[2].push({t: date, y: Math.round(sol_history[i][3])}); //I0
        data_cat[3].push({t: date, y: Math.round((1-c)*sol_history[i][4])}); //mild hidden: (1-c)*I1
        data_cat[4].push({t: date, y: Math.round(sol_history[i][9])}); //R hidden
        data_cat[5].push({t: date, y: Math.round(sol_history[i][8])}); //R diagnosed
        data_cat[6].push({t: date, y: Math.round(c*sol_history[i][4]) + Math.round(sol_history[i][5])}); //mild diagnosed: c*I1 + Iq
        data_cat[7].push({t: date, y: Math.round(sol_history[i][6])}); //I2
        data_cat[8].push({t: date, y: Math.round(sol_history[i][7])}); //I3
        data_cat[9].push({t: date, y: Math.round(sol_history[i][10])}); //D

        let num_confirmed_cases = 0;
        for (let j = 0; j < report_sum_indices.length; ++j)
          num_confirmed_cases += Math.round(sol_history[i][report_sum_indices[j]]);
        num_confirmed_cases += Math.round(c*sol_history[i][4]); //c*I1
        data_agg.push({t: date, y: num_confirmed_cases});
      }

      return {aggregated: data_agg, categorized: data_cat};
    }

    function predictModel(params)
    {
      //Periods [days]
      const T_incub0 = 3;
      const T_incub1 = 2;
      const T_asympt = 6;
      const T_mild   = 6;
      const T_severe = 4;
      const T_icu    = 10;

      //Probabilities
      // const prob_E0_E1 = 1;  //non-infectious exposed to infectious exposed
      // const prob_E1_I0 = 0.3 / prob_E0_E1;  //exposed to asymptomatic
      // const prob_E1_I1 = 1 - prob_E1_I0;  //exposed to mild
      // const prob_I0_R  = 1;
      // const prob_I1_R   = 0.51/(prob_E0_E1*prob_E1_I1); //mild to recovered
      // const prob_I1_I2  = 1 - prob_I1_R; //mild to severe
      // const prob_I2_R   = 0.14/(prob_E0_E1*prob_E1_I1*prob_I1_I2); //severe to recovered
      // const prob_I2_I3  = 1 - prob_I2_R; //severe to critical
      // const prob_I3_D   = 0.02/(prob_E0_E1*prob_E1_I1*prob_I1_I2*prob_I2_I3); //critical to dead
      // const prob_I3_R   = 1 - prob_I3_D; //critical to recovered

      const prob_E0_E1 = 1;  //non-infectious exposed to infectious exposed
      const prob_E1_I0 = 0.3 / prob_E0_E1;  //exposed to asymptomatic
      const prob_E1_I1 = 1 - prob_E1_I0;  //exposed to mild
      const prob_I0_R  = 1;
      const prob_I1_R   = 0.8; //mild to recovered
      const prob_I1_I2  = 1 - prob_I1_R; //mild to severe
      const prob_I2_R   = 0.15/(prob_I1_I2); //severe to recovered
      const prob_I2_I3  = 1 - prob_I2_R; //severe to critical
      const prob_I3_D   = 0.02/(prob_I1_I2*prob_I2_I3); //critical to dead
      const prob_I3_R   = 1 - prob_I3_D; //critical to recovered


      //Rates [1/day]
      const a0 = (1/T_incub0) * prob_E0_E1;
      const a10 = (1/T_incub1) * prob_E1_I0;
      const a11 = (1/T_incub1) * prob_E1_I1;
      const a1 = a10 + a11;
      const g0 = (1/T_asympt) * prob_I0_R;
      const g1 = (1/T_mild)   * prob_I1_R;
      const p1 = (1/T_mild)   * prob_I1_I2;
      const g2 = (1/T_severe) * prob_I2_R;
      const p2 = (1/T_severe) * prob_I2_I3;
      const g3 = (1/T_icu)    * prob_I3_R;
      const mu = (1/T_icu)    * prob_I3_D;

      let N = 21.4e6; //Population of Sri Lanka

      //Initial solution
      let E0_0 = params.E0_0;
      let E1_0 = 0;
      let I0_0 = 0;
      let I1_0 = 0;
      let Iq_0 = 0;
      let I2_0 = 0;
      let I3_0 = 0;
      let Rd_0 = 1;  //One patient had already recovered in SL
      let Rh_0 = 0;
      let D_0 = 0;
      let S_0 = N - E0_0 - E1_0 - I0_0 - I1_0 - Iq_0 - I2_0 - I3_0 - Rd_0 - Rh_0 - D_0;

      //Solution vector: [S, E0, E1, I0, I1, Iq, I2, I3, Rd, Rh, D]
      let solution_hist = [[S_0, E0_0, E1_0, I0_0, I1_0, Iq_0, I2_0, I3_0, Rd_0, Rh_0, D_0]];

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
          let dS = -(b1*(u[2] + u[3] + u[4]) + b2*u[6] + b3*u[7])*u[0];
          let du = [ dS,                                                                                          //S
                    -dS - a0*u[1],                                                                                //E0
                          a0*u[1] - a1*u[2],                                                                      //E1
                                   a10*u[2] - g0*u[3],                                                            //I0
                                   a11*u[2]           - (g1 + p1)*u[4],                                           //I1
                      q_input                                - (g1 + p1)*u[5],                                    //Iq
                                                              p1*(u[4] + u[5]) - (g2 + p2)*u[6],                  //I2
                                                                                        p2*u[6] - (g3 + mu)*u[7], //I3
                                                            g1*(c*u[4] + u[5])        + g2*u[6]        + g3*u[7], //Rd
                                              g0*u[3] +  (1-c)*g1*u[4],                                           //Rh
                                                                                                         mu*u[7]  //D
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


    function optimizeParameters()
    {
      let params = initializeSimulationParameters(data_SL.length, 0); //no prediction

      let res = getFitResidual(params);
      let resnorm = getL2Norm(res);
      let resnorm_init = resnorm;

      let m = res.length;
      let n_b1 = (params.T_hist - 1); //no. of beta_1 values to optimize
      let n_c = (params.T_hist - 1); //no. of c values to optimize
      let n = n_b1 + 0*n_c; //beta_1 values, E0_0, I1h_0

      let param_vec = new Array(n).fill(0);
      for (let i = 0; i < n_b1; ++i)
        param_vec[i] = params.b1N[i];
      // for (let i = 0; i < n_c; ++i)
      //   param_vec[n_b1 + i] = params.diag_frac[i];

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
            else
            {
              // if (param_vec[i] < 0.0)
              //   param_vec[i] = 0.0;
              // else if (param_vec[i] > 1.0)
              //   param_vec[i] = 1.0;
              // params.diag_frac[i-n_b1] = param_vec[i];
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

      console.log("Final resnorm: " + resnorm + ", rel: " + resnorm/resnorm_init);

      for (let i = 0; i < n_b1; ++i)
        sim_params.b1N[i] = params.b1N[i];

      // for (let i = 0; i < n_c; ++i)
      //     sim_params.diag_frac[i] = params.diag_frac[i];

      //Extend last beta and c value to future
      for (let i = n_b1; i < sim_params.b1N.length; ++i)
        sim_params.b1N[i] = sim_params.b1N[n_b1-1];

      // for (let i = n_c; i < sim_params.diag_frac.length; ++i)
      //   sim_params.diag_frac[i] = sim_params.diag_frac[n_c-1];

      updateParameters();
      return params;
    }

    function getFitResidual(params)
    {
      const num_eq = 2;
      let sol_hist = predictModel(params);
      let residual = new Array(num_eq*(sol_hist.length-1)).fill(0);

      for (let i = 1; i < sol_hist.length; ++i)
      {
        let weight0 = (data_raw_SL[i].y[0] == 0) ? 1 : (1/data_raw_SL[i].y[0]);
        let weight1 = (data_raw_SL[i].y[1] == 0) ? 1 : (1/data_raw_SL[i].y[1]);
        let weight2 = (data_raw_SL[i].y[2] == 0) ? 1 : (1/data_raw_SL[i].y[2]);

        let I_diag = params.diag_frac[i]*sol_hist[i][4] + sol_hist[i][5] + sol_hist[i][6] + sol_hist[i][7];

        //I1d + I1q + I2 + I3 = no. of diagnosed patients in hospitals
        //residual[num_eq*(i-1)] = (I_diag - data_raw_SL[i].y[0]);
        residual[num_eq*(i-1)] = (I_diag + sol_hist[i][8] - data_raw_SL[i].y[0] - data_raw_SL[i].y[1]);
        //residual[num_eq*(i-1)] = I_diag + sol_hist[i][8] + sol_hist[i][10] - data_raw_SL[i].y[0] - data_raw_SL[i].y[1] - data_raw_SL[i].y[2];

        //Rd = no. of recovered patients
        //residual[num_eq*(i-1) + 1] = (sol_hist[i][8] - data_raw_SL[i].y[1]) * weight1;

        //D = no. of fatalities
        residual[num_eq*(i-1) + 1] = (sol_hist[i][10] - data_raw_SL[i].y[2]) * 10;
      }
      return residual;
    }

    function getFitJacobian(params)
    {
      let m = 2*(params.T_hist - 1);
      let n_b1 = (params.T_hist - 1); //no. of beta_1 values to optimize
      let n_c = (params.T_hist - 1); //no. of c values to optimize
      let n = n_b1 + 0*n_c; //beta_1 values, E0_0, I1h_0
      let jac = new Array(m*n).fill(0);

      const delta_b1 = 1e-5;
      const delta_c = 1e-5;

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

      // for (let j = 0; j < n_c; ++j)
      // {
      //   //Compute finite difference
      //   params.diag_frac[j] += delta_c;
      //   let Rp = getFitResidual(params);
      //   params.diag_frac[j] -= 2*delta_c;
      //   let Rm = getFitResidual(params);
      //   params.diag_frac[j] += delta_c;
      //
      //   for (let i = 0; i < m; ++i)
      //     jac[i*n + n_b1 + j] = (Rp[i] - Rm[i])/(2*delta_c);
      // }

      return jac;
    }

    function getL2Norm(vec)
    {
      let norm = 0.0;
      for (let i = 0; i < vec.length; ++i)
        norm += vec[i]*vec[i];
      return Math.sqrt(norm);
    }

});