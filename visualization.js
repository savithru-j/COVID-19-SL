// COVID-19 Simulation Tool in JavaScript
// Copyright 2020 Savithru Jayasinghe and Dushan Wadduwage
// Licensed under the MIT License (LICENSE.txt)

'use strict';

var date_format = 'YYYY-MM-DD';

var country_data;
var district_data;

//Chart objects
var country_chart;
var district_chart;

var countries_to_plot = ['Sri Lanka', 'United States', 'United Kingdom', 'India',
                         'Pakistan','China', 'Singapore', 'Japan', 'Brazil', 'Spain', 'Italy',
                         'Oman','Thailand'];

window.onload = function()
{
  loadCountryData();
  setupCountryChart();
  setupDistrictChart();
}

function loadCountryData()
{
  let case_threshold = 100;

  country_data = [];
  for (let country_name of Object.keys(world_data))
  {
      let data_array = world_data[country_name];
      let processed_data = [];

      let start_date = null;
      for (let data of data_array)
      {
        if (data.confirmed >= case_threshold)
        {
          if (start_date == null)
            start_date = moment(data.date, date_format);

          let current_date = moment(data.date, date_format);
          let t_elapsed = current_date.diff(start_date, 'days');
          processed_data.push({x:t_elapsed, y: data.confirmed});
        }
      }
      if (country_name == 'US')
        country_name = 'United States';
      country_data[country_name] = processed_data;
  }
}

function setupCountryChart()
{
  let chart_element = document.getElementById('chart_country');

  let dataseries_vec = [];
  let xmax = 0;
  for (let country_name of countries_to_plot)
  {
    let data_hist = country_data[country_name];
    dataseries_vec.push({
      name: country_name,
      data: data_hist
    });
    xmax = Math.max(xmax, data_hist[data_hist.length-1].x);
  }

  var options = {
    chart: {
      type: 'line',
      animations: {
        enabled: true,
        easing: 'easeinout',
        speed: 500,
        animateGradually: {
            enabled: false,
        },
        dynamicAnimation: {
            enabled: true,
            speed: 350
        }
      }
    },
    series: dataseries_vec,
    xaxis: {
      type: 'numeric',
      title: {
        text: 'No. of days elapsed since 100 cases',
        style: {
            fontSize: "16px",
            color: '#666'
        },
        offsetY: 10
      },
    },
    yaxis: {
      logarithmic: true,
      title: {
        text: 'No. of cases',
        style: {
            fontSize: "16px",
            color: '#666',
        },
      },
    },
    legend: {
      offsetY: 15,
      position: 'right'
    },
    grid: {
      xaxis: {
        lines: {
          show: true
        }
      },
      yaxis: {
        lines: {
          show: true
        }
      }
    },
    colors: ['#1f78b4','#a6cee3','#b2df8a','#33a02c','#fb9a99','#e31a1c',
             '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928','#ffff99',
             '#9665c7'],
    fill: {
      type: 'solid'
    },
    stroke: {
      width: 1.5,
      curve: 'smooth',
    },
  }

  country_chart = new ApexCharts(chart_element, options);
  country_chart.render();
};

function setupDistrictChart()
{
  let chart_element = document.getElementById('chart_district');

  let dataseries_vec = [];
  // let districts = ['Colombo','Gampaha','Jaffna'];
  for (let district_name of Object.keys(SL_district_data))
  {
    dataseries_vec.push({
      name: district_name,
      data: SL_district_data[district_name]
    });
  }

  let yaxis_log = document.querySelector('input[name="plot_yaxis_type"]:checked').value == 'logarithmic';

  var options = {
    chart: {
      type: 'line',
      animations: {
        enabled: true,
        easing: 'easeinout',
        speed: 500,
        animateGradually: {
            enabled: false,
        },
        dynamicAnimation: {
            enabled: true,
            speed: 350
        }
      }
    },
    series: dataseries_vec,
    xaxis: {
      type: 'datetime',
      title: {
        text: 'Date',
        style: {
            fontSize: "16px",
            color: '#666'
        },
        offsetY: 10
      },
    },
    yaxis: {
      logarithmic: yaxis_log,
      title: {
        text: 'No. of cases',
        style: {
            fontSize: "16px",
            color: '#666',
        },
      },
    },
    legend: {
      offsetY: 15,
      position: 'right'
    },
    grid: {
      xaxis: {
        lines: {
          show: true
        }
      },
      yaxis: {
        lines: {
          show: true
        }
      }
    },
    colors: ['#1f78b4','#a6cee3','#b2df8a','#33a02c','#fb9a99','#e31a1c',
             '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928','#ffff99',
             '#9665c7'],
    fill: {
      type: 'solid'
    },
    stroke: {
      width: 1.5,
      curve: 'smooth',
    },
    tooltip: {
      enabled: true,
      shared: false
    }
  }

  district_chart = new ApexCharts(chart_element, options);
  district_chart.render();
};

function updateYAxisType()
{
  let yaxis_log = document.querySelector('input[name="plot_yaxis_type"]:checked').value == 'logarithmic';
  district_chart.updateOptions({
    yaxis: {
      logarithmic: yaxis_log,
    }
  });
}
