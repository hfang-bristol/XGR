(function() {
	function barChart (element, model) {
		this.element = element;
		this.model = model;
	}

	barChart.prototype.render = function() {
		
		// options
		var options = {
			chart: {
				renderTo: this.element,
				type: 'bar',
				reflow: false,
			},
			title: {
				text: '',
			},
			subtitle: {
				text: '',
			},
			xAxis: {
				categories: [],
				lineWidth: 2,
				lineColor: '#000000',
				labels: {
					overflow: 'justify',
					style: {
						color: 'black',
						fontSize:'13px'
					}
				}
			},
			yAxis: {
				opposite: true,
				lineWidth: 2,
				lineColor: '#000000',
				gridLineWidth: 1,
				plotLines: [
					{
						color: "rgba(0,100,0,0.8)",
						width: 2,
						value: 1.3
            		},
            	],
				min: 0,
				//tickInterval: 2,
				title: {
					text: '-log10(FDR)',
					align: 'middle',
					style: {
						color: 'black',
						fontSize:'20px'
					}
				},
				labels: {
					overflow: 'justify',
					style: {
						color: 'black',
						fontSize:'15px'
					}
				}
			},
	
			tooltip: {
				useHTML: true,
                formatter: function() {
                	var fdr=Math.pow(10,-1*this.y).toPrecision(2);
                	var Zscore=this.point.Zscore;
                	var Pvalue=this.point.Pvalue;
                	
                	var nAnno=this.point.nAnno;
                	var nOverlap=this.point.nOverlap;
                	var Member=this.point.Member;
                	
                	return '<b>'+ this.x + '</b><br/>' +
                	'<li><span style="color:'+ this.point.series.color +'">FDR</span> = ' + fdr +
                	'<li><span style="color:'+ this.point.series.color +'">P-value</span> = ' + Pvalue +
                	'<li><span style="color:'+ this.point.series.color +'">Z-score</span> = ' + Zscore +
                	'<li><span style="color:'+ this.point.series.color +'">nOverlaps/nAnnos</span> = ' + nOverlap + '/' + nAnno +
                	'<li><span style="color:'+ this.point.series.color +'">Members</span>:<br\>' + 
                	'<span style="font-size:8px; display:inline-block; max-width:300px; white-space:normal; overflow:hidden; text-overflow:ellipsis;">' + Member + '</span>'
                	;
                	
                },
				crosshairs: {
					color: 'gray'
				},
				//backgroundColor: '#FCFFC5',
				backgroundColor: 'rgba(252, 255, 197, 0.9)',
				borderColor: 'gray',
				borderRadius: 10,
				borderWidth: 3,
				
			},
			plotOptions: {
				bar: {
					dataLabels: {
						enabled: false
					}
				}
			},
			legend: {
				layout: 'vertical',
				align: 'right',
				verticalAlign: 'bottom',
				x: -40,
				y: 80,
				floating: true,
				borderWidth: 1,
				backgroundColor: '#FFFFFF',
				shadow: true
			},
			credits: {
				enabled: false
			},
			series: []
		}
		
		// getJSON
        $.getJSON(this.model, function(json) {
            
            	var data=json['data'];
            	
            	// sort according to adjp
            	data.sort(function(a, b) {
					return a['adjp'] - b['adjp'];
				});
            	
            	var myCategories = [];
            	var myData = [];
				
				// replace 0 with nonzero min
				var nonzero_min = 1;
				for(i=0; i<data.length; ++i){
					if(data[i]['adjp']!=0){
						nonzero_min=data[i]['adjp'];
						break;
					}
				}
				
           		var i;
           		var i_max=30;
            	for(i=0; i<data.length; ++i){
            		if(data[i]['adjp']==0){
            			data[i]['adjp']=nonzero_min;
            		}
            		var tmpY=-1*Math.log10(data[i]['adjp']);
            		
            		var tmpName=data[i]['name'];
            		var tmpZ=data[i]['zscore'];
            		var tmpPvalue=data[i]['pvalue'];
            		var tmpFDR=data[i]['adjp'];
            		
            		var tmpAnno=data[i]['nAnno'];
            		var tmpOverlap=data[i]['nOverlap'];
            		
            		var tmpMember=data[i]['members'];
            		
            		//if(i<i_max && tmpY > (-1*Math.log10(0.05))){
            		if(i<i_max){
						myCategories.push(tmpName);
						
						myData.push({
							y: tmpY,
							Zscore: tmpZ,
							Pvalue: tmpPvalue,
							FDR: tmpFDR,
							nAnno: tmpAnno,
							nOverlap: tmpOverlap,
							Member: tmpMember,
						});
					}
            	}
            	
            	var mySeries = {
            		name: "Terms",
    				data: myData
            	}
    
                options.xAxis.categories = myCategories;
                options.series[0] = mySeries;
                var chart = new Highcharts.Chart(options);
        });
	
	}

	window.barChart = barChart;
})();