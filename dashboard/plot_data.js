
'use strict';

var date_format = 'DD-MM-YYYY';

function getDataSriLanka()
{
  let data = [];
  data.push({t: moment('10-03-2020', date_format), y: 1});
  data.push({t: moment('11-03-2020', date_format), y: 2});
  data.push({t: moment('12-03-2020', date_format), y: 3});
  data.push({t: moment('13-03-2020', date_format), y: 6});
  data.push({t: moment('14-03-2020', date_format), y: 11});
  data.push({t: moment('15-03-2020', date_format), y: 19});
  data.push({t: moment('16-03-2020', date_format), y: 29});
  data.push({t: moment('17-03-2020', date_format), y: 42});
  data.push({t: moment('18-03-2020', date_format), y: 53});
  data.push({t: moment('19-03-2020', date_format), y: 66});
  data.push({t: moment('20-03-2020', date_format), y: 72});
  data.push({t: moment('21-03-2020', date_format), y: 78});
  data.push({t: moment('22-03-2020', date_format), y: 87});
  data.push({t: moment('23-03-2020', date_format), y: 97});
  data.push({t: moment('24-03-2020', date_format), y: 102});
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

function getDataStacked()
{
  let data0 = [];
  data0.push({t: moment('10-03-2020', date_format), y: 20});
  data0.push({t: moment('11-03-2020', date_format), y: 30});
  data0.push({t: moment('12-03-2020', date_format), y: 40});
  
  let data1 = [];
  data1.push({t: moment('10-03-2020', date_format), y: 40});
  data1.push({t: moment('11-03-2020', date_format), y: 60});
  data1.push({t: moment('12-03-2020', date_format), y: 80});
  
  let data2 = [];
  data2.push({t: moment('10-03-2020', date_format), y: 60});
  data2.push({t: moment('11-03-2020', date_format), y: 50});
  data2.push({t: moment('12-03-2020', date_format), y: 90});
  return [data0, data1, data2];
}

var dataSL = getDataSriLanka();
var dataIT = getDataItaly();

var sim_params = initializeSimulationParameters(dataSL.length);

var dataSLpredstack = [];
var dataSLpredicted = getPredictionData(dataSL[0].t);

var chart = [];
var chart_config = [];

window.onload = function()
{
  updateChart();
}

function updateChart()
{
  let canvas = document.getElementById('chart_canvas');
  canvas.width = 0.7*window.innerWidth;
  canvas.height = 0.65*window.innerHeight;
  let ctx = canvas.getContext('2d');
    
  let div_controls = document.getElementById('chart_controls');
  div_controls.style.width = (0.15*window.innerWidth) + "px";
  
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
  
	chart_config = {
			data: {
				datasets: [{
					label: 'Sri Lanka - confirmed cases',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#66b3ff',
					data: dataSL,
					type: 'line',
					fill: true,
					borderWidth: 2
				},
				{
					label: 'Sri Lanka - predicted',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#66b3ff',
					borderDash: [5, 5],
					data: dataSLpredicted,
					type: 'line',
					fill: true,
					borderWidth: 2
				},
				{
					label: 'Italy - confirmed cases',
					backgroundColor: 'rgba(1,1,1,0)',
					borderColor: '#fcaf3e',
					data: dataIT,
					type: 'line',
					fill: true,
					borderWidth: 2,
					hidden: true
				},
				{
					label: 'Susceptible',
					backgroundColor: '#d3d7cf',
					data: dataSLpredstack[0],
					type: 'bar',
					stack: 'stack0',
					hidden: true
				},
				{
					label: 'Exposed',
					backgroundColor: '#e9b96e',
					data: dataSLpredstack[1],
					type: 'bar',
					stack: 'stack0'
				},
				{
					label: 'Mildly infected',
					backgroundColor: '#edd400',
					data: dataSLpredstack[2],
					type: 'bar',
					stack: 'stack0'
				},
				{
					label: 'Severely infected',
					backgroundColor: '#f57900',
					data: dataSLpredstack[3],
					type: 'bar',
					stack: 'stack0'
				},
				{
					label: 'Critically infected',
					backgroundColor: '#cc0000',
					data: dataSLpredstack[4],
					type: 'bar',
					stack: 'stack0'
				},
				{
					label: 'Recovered',
					backgroundColor: '#73d216',
					data: dataSLpredstack[5],
					type: 'bar',
					stack: 'stack0'
				},
				{
					label: 'Fatal',
					backgroundColor: '#555753',
					data: dataSLpredstack[6],
					type: 'bar',
					stack: 'stack0'
				}
				]
			},
			options: {
			  responsive: false,
        maintainAspectRatio: false,
				scales: {
					xAxes: [xaxis_config],
					yAxes: [yaxis_config]
				},
				legend: {
            display: true,
            boxWidth: 10,
        },
//        tooltips: {
//          callbacks: {
//                label: function(tooltipItem, data) {
//                    //var type = data.datasets[tooltipItem.datasetIndex].label;
//                    //var value = data.datasets[tooltipItem.datasetIndex].data[tooltipItem.index].y;

//                    // Loop through all datasets to get the actual total of the index
//                    var total = 0;
//                    let labels = [];
//                    for (var i = 3; i < data.datasets.length; i++)
//                    {
//                      labels.push(data.datasets[i].label + ": " + data.datasets[i].data[tooltipItem.index].y);
//                      total += data.datasets[i].data[tooltipItem.index].y;
//                    }
//                    labels.push("Total: " + total);
//                    return labels;
//                }
//            }
//        },
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
    max: 100000,
    callback: function (value, index, values) {
      if (value === 1000000) return "1000000";
      if (value === 100000) return "100000";
      if (value === 10000) return "10000";
      if (value === 1000) return "1000";
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
  chart_config.data.datasets[2].hidden = !overlay;
  chart.update();
}

function alignTimelines(align)
{
  if (align)
  {
    for (let datapair of chart_config.data.datasets[2].data)
      datapair.t.add(41,'days'); //Add offset to Italy dataset
  }
  else //load original datasets
  {
    chart_config.data.datasets[2].data = getDataItaly();
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
  const total_length = hist_length + pred_length;
  let params = {
    T_hist: hist_length,
    T_pred: pred_length,
    dt: 0.5/24.0,        //timestep size [days]
    b1N: new Array(total_length).fill(0.5), //transmission rate from mild to susceptible
    b2N: new Array(total_length).fill(0.0), //transmission rate from severe to susceptible
    b3N: new Array(total_length).fill(0.0)  //transmission rate from critical to susceptible
  }
  return params;
}

function updateParameters()
{
  let requires_update = false;
  
  let slider_element_ids = ["slider_b1", "slider_b2", "slider_b3"];
  let param_arrays = [sim_params.b1N, sim_params.b2N, sim_params.b3N];
  
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
        }
    }
  }
  
//  let slider_dt = document.getElementById("slider_dt");
//  if (slider_dt)
//  {
//    let val = 1.0/Math.pow(2,Number(slider_dt.value));
//    if (sim_params.dt != val)
//    {
//      sim_params.dt = val;
//      requires_update = true;
//    }
//  }
  
  if (requires_update)
  {
    dataSLpredicted = getPredictionData(dataSL[0].t);
    chart_config.data.datasets[1].data = dataSLpredicted;
    chart.update();
  }
}

function getPredictionData(start_date)
{
  let data = [];
  let sol_history = predictModel(sim_params);
  
  dataSLpredstack = new Array(7).fill([]);
  
  for (let i = 0; i < sol_history.length; i++)
  {
    let date = start_date.clone().add(i,'days');
    
    let num_cases = 0;
    for (let j = 1; j < 7; ++j)
      num_cases += sol_history[i][j];
      
    data.push({t: date, y: Math.round(num_cases)});
    
    //for (let j = 0; j < 7; ++j)
    //  dataSLpredstack[j].push({t: date, y: Math.round(sol_history[i][j])});
  }
  
  return data;
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
  let a   = (1/T_incub)  *prob_I1_E;
  let g1  = (1/T_mild)   * prob_R_I1;
  let p1  = (1/T_mild)   * prob_I2_I1;
  let g2  = (1/T_severe) * prob_R_I2;
  let p2  = (1/T_severe) * prob_I3_I2;
  let g3  = (1/T_icu)    * prob_R_I3;
  let u   = (1/T_icu)    * prob_D_I3;

  //Transmission rates: beta values are always scaled by the population N

  let b1N = 0.5; //rate at which mildly infected individuals transmit to susceptible people
  let b2N = 0;   //rate at which severely infected individuals transmit to susceptible people
  let b3N = 0;   //rate at which critically infected individuals transmit to susceptible people
  
  let N = 21.4e6; //Population of Sri Lanka
   
  //Initial solution
  let E0 = 20;
  let S0 = N - E0;
  let I1_0 = 0;
  let I2_0 = 0;
  let I3_0 = 0;
  let R0 = 0;
  let D0 = 0;
  
  //Solution vector: [S, E, I1, I2, I3, R, D] 
  let solution_hist = [[S0, E0, I1_0, I2_0, I3_0, R0, D0]];
  
  const nt = params.T_hist + params.T_pred;
  const nt_sub = 1.0/params.dt;
  
  for (let i = 0; i < nt; i++) 
  {
    let sol = solution_hist[i].slice(); //copy last solution vector
    
    let b1 = params.b1N[i] / N;
    let b2 = params.b2N[i] / N;
    let b3 = params.b3N[i] / N;
    
    for (let j = 0; j < nt_sub; j++)
    {
      let dS = -(b1*sol[2] + b2*sol[3] + b3*sol[4])*sol[0]; 
      let dsol = [dS,                                //dS
                  -dS - a*sol[1],                    //dE
                  a*sol[1] - (g1 + p1)*sol[2],       //dI1
                  p1*sol[2] - (g2 + p2)*sol[3],      //dI2
                  p2*sol[3] - (g3 + u)*sol[4],       //dI3
                  g1*sol[2] + g2*sol[3] + g3*sol[4], //dR
                  u*sol[4]                           //dD
                 ];
      
      for (let k = 0; k < sol.length; k++)
        sol[k] += dsol[k]*params.dt;
    } //sub-timestepping [hrs]

    solution_hist.push(sol); //save solution daily
  }
  
  return solution_hist;
}
