
'use strict';

var date_format = 'DD-MM-YYYY';

function getDataSriLanka()
{
  var data = [];
  data.push({t: moment('11-03-2020', date_format), y: 1});
  data.push({t: moment('12-03-2020', date_format), y: 2});
  data.push({t: moment('13-03-2020', date_format), y: 3});
  data.push({t: moment('14-03-2020', date_format), y: 6});
  data.push({t: moment('15-03-2020', date_format), y: 11});
  data.push({t: moment('16-03-2020', date_format), y: 19});
  data.push({t: moment('17-03-2020', date_format), y: 29});
  data.push({t: moment('18-03-2020', date_format), y: 42});
  data.push({t: moment('19-03-2020', date_format), y: 53});
  data.push({t: moment('20-03-2020', date_format), y: 66});
  data.push({t: moment('21-03-2020', date_format), y: 72});
  data.push({t: moment('22-03-2020', date_format), y: 78});
  return data;
}

function getDataItaly()
{
  var data = [];
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
  return data;
}

var chart = [];
var chart_config = [];

window.onload = function()
{
  updateChart();
}


function updateChart()
{
  var canvas = document.getElementById('chart_canvas');
  canvas.width = 0.7*window.innerWidth;
  canvas.height = 0.65*window.innerHeight;
  var ctx = canvas.getContext('2d');
    
  var div_controls = document.getElementById('chart_controls');
  div_controls.style.width = (0.15*window.innerWidth) + "px";
  
  var check_logy = document.getElementById('check_log_y');
 
  var dataSL = getDataSriLanka();
  var dataIT = getDataItaly();
  
  var xaxis_config = {
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
			}
		};
		
  var yaxis_config = {
			gridLines: {
				drawBorder: false
			},
			type: 'linear',
			scaleLabel: {
				display: true,
				labelString: 'No. of confirmed cases',
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
					label: 'Sri Lanka',
					borderColor: '#66b3ff',
					data: dataSL,
					type: 'line',
					fill: false,
					borderWidth: 2
				},
				{
					label: 'Italy',
					borderColor: '#fcaf3e',
					data: dataIT,
					type: 'line',
					fill: false,
					borderWidth: 2,
					hidden: true
				}]
			},
			options: {
			  responsive: false,
        maintainAspectRatio: false,
				scales: {
					xAxes: [xaxis_config],
					yAxes: [yaxis_config]
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
  chart_config.data.datasets[1].hidden = !overlay;
  chart.update();
}

function alignTimelines(align)
{
  if (align)
  {
    for (let datapair of chart_config.data.datasets[1].data)
      datapair.t.add(41,'days'); //Add offset to Italy dataset
  }
  else //load original datasets
  {
    chart_config.data.datasets[1].data = getDataItaly();
  }
  chart.update();
}

function resetZoom()
{
  chart.resetZoom();
}
