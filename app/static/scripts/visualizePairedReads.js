function visualizePairedReads(json_paired_data, div = 'pairsImage') {
	var pairsChart = new CanvasJS.Chart(div , {
		exportEnabled: true,
		animationEnabled: true,
		theme: "light2",
		title:{
			text: "Percentage Paired reads"
		},
		legend:{
			cursor: "pointer"
		},
		data: [{
			type: "pie",
			showInLegend: true,
			toolTipContent: "{name}: <strong>{y}%</strong>",
			indexLabel: "{name} - {y}%",
			dataPoints: json_paired_data
		}]
	});
	pairsChart.render();
	enable_all_buttons();
};
