function visualizeNucleotideAndQualityPercentages(json_data) {
    var nucleotidePercentagesChart = new CanvasJS.Chart("sequenceTrimImage", {
        exportEnabled: true,
        animationEnabled: true,
        zoomEnabled: true,
        theme: "light2",
        title: {
            text: "Nucleotide percentage per position"
        },
		axisY: {
			labelFormatter: function(e) {
				return Math.round((e.value / json_data["max"]) * 100) + " %";
            },
            minimum: 0,
            maximum: json_data["max"]
		},
        legend: {
            cursor: "pointer",
        },
        data: [{
            name: "A",
            type: "scatter",
            showInLegend: true,
            yValueFormatString: "#0.## reads",
            dataPoints: json_data["A"]
        },
        {
            name: "C",
            type: "scatter",
            showInLegend: true,
            yValueFormatString: "#0.## reads",
            dataPoints: json_data["C"]
        },
        {
            name: "G",
            type: "scatter",
            showInLegend: true,
            yValueFormatString: "#0.## reads",
            dataPoints: json_data["G"]
        },
        {
            name: "T",
            type: "scatter",
            showInLegend: true,
            yValueFormatString: "#0.## reads",
            dataPoints: json_data["T"]
        }]
    });

    var sequenceQualityChart = new CanvasJS.Chart("sequenceQualityImage", {
        exportEnabled: true,
        animationEnabled: true,
        zoomEnabled: true,
        theme: "light2",
        title:{
            text: "Sequence quality per position"
        },
        legend: {
            cursor: "pointer"
        },
        data: [{
            type: "rangeArea",
            showInLegend: true,
            name: "Range",
            color: "#0080FF",
            dataPoints: json_data["range"]
        },
        {
            type: "line",
            showInLegend: true,
            name: "Average",
            color: "#FF4040",
            dataPoints: json_data["average"]
        }]
    });

    nucleotidePercentagesChart.render();
    sequenceQualityChart.render();
    enable_all_buttons();
};
