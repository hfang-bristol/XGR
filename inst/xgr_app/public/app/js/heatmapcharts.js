(function() {
	function heatmapChart (element, model) {
		this.element = element;
		this.model = model;
	}

	heatmapChart.prototype.render = function() {
		
		// options
		var options = {
			chart: {
				renderTo: this.element,
				type: 'heatmap',
				reflow: false,
				marginTop: 40,
				marginBottom: 80,
				plotBorderWidth: 1
			},
			title: {
				text: '',
			},
			subtitle: {
				text: '',
			},
			xAxis: {
				categories: [],
				labels: {
					overflow: 'justify',
					style: {
						color: Highcharts.getOptions().colors[0],
						fontSize:'12px'
					}
				}
			},
			yAxis: {
				categories: [],
				labels: {
					overflow: 'justify',
					style: {
						color: Highcharts.getOptions().colors[0],
						fontSize:'12px'
					}
				},
				title: null
			},
			
			colorAxis: {
				min: 0,
				minColor: '#FFFFFF',
				maxColor: '#FF0000',
				//maxColor: Highcharts.getOptions().colors[0],
			},
			legend: {
				align: 'right',
				layout: 'vertical',
				margin: 0,
				verticalAlign: 'top',
				y: 25,
				//symbolHeight: 280
			},
			
			plotOptions: {
					series: {
						turboThreshold:0
					}
			},
			
			tooltip: {
				formatter: function () {
					return '<b>' + this.series.xAxis.categories[this.point.x] + '</b> and <b>' + this.series.yAxis.categories[this.point.y] + '</b>:<br>Similarity = ' + this.point.value;
				}
			},
			credits: {
				enabled: false
			},
			series: []
		}
		
		// getJSON
        $.getJSON(this.model, function(json) {
            
            	var data=json['data'];
            	var myCategories=json['category'];
				
				var mySeries=[{
						borderWidth: 1,
						data: data,
						dataLabels: {
							enabled: false,
							color: '#000000'
						}
					}];
					
                options.xAxis.categories = myCategories;
                options.yAxis.categories = myCategories;
                options.series = mySeries;
                var chart = new Highcharts.Chart(options);
                
	    // Horrible hack necessary to display the chart properly
		setTimeout(function () {
			$(window).trigger('resize');
		}, 0);
                
                
$(window).resize(function() {
    width = this.element.width();
    height = width;
    chart.setSize(width, height, doAnimation = true);
});
                
        });
	
	}

	window.heatmapChart = heatmapChart;
	
})();