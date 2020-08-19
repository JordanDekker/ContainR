function visualizeSequenceLength(data) {
    var sequenceChart = new CanvasJS.Chart("sequenceImage", {
        animationEnabled: true,
        exportEnabled: true,
        zoomEnabled: true,
        theme: "light2",
        title: {
            text: "Distribution of sequence length"
        },
        axisX: {
            title: "Sequence length (nucleotides)",
            labelFormatter: function(){
                return " ";
              }
        },
        axisY: {
            title: "Count (10 log)",
            logarithmic:  true
        },
        legend: {
            cursor: "pointer",
            verticalAlign: "bottom",
            horizontalAlign: "left",
            dockInsidePlotArea: false
        },
        data: [{
            type: "column",
            showInLegend: true,
            name: "Forward reads",
            indexLabel: "{x}",
		    indexLabelFontColor: "#000000",
		    indexLabelPlacement: "outside",
            dataPoints: data["fw_seq_length"],
            toolTipContent: "{y} reads with {x} nucleotides"
        },
            {
                type: "column",
                showInLegend: true,
                name: "Reverse reads",
                indexLabel: "{x}",
		        indexLabelFontColor: "#000000",
		        indexLabelPlacement: "outside",
                dataPoints: data["rv_seq_length"],
                toolTipContent: "{y} reads with {x} nucleotides"
            }]
    });
    sequenceChart.render();
    enable_all_buttons();
}