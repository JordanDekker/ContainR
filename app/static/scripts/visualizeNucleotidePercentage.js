function visualizeNucleotidePercentages(nucleotidePercentagesData, div, title="Clustered Nucleotide Percentage") {
    // TODO Add text to boxplots explaining what the whiskers mean: EG min and max values still within IQR
    outlierdata = []
    for (const i of nucleotidePercentagesData) {
        for (const outlier of i["outliers"]) {
            outlierdata.push({"x": i["x"], "y": outlier})
        }
    }

    var nucleotideChart = new CanvasJS.Chart(div, {
        animationEnabled: true,
        theme: "light2", // "light1", "light2", "dark1", "dark2"
        title:{
            text: title
        },
        subtitles: [{
            text: "",
            fontSize: 15
        }],
        axisY: {
            title: "Percentage of nucleotide in reads",
            includeZero: false,
            tickLength: 0,
            gridDashType: "dash",
            stripLines: [{
                value: 25,
                label: "25%",
                labelFontColor: "#70e5ea",
                showOnTop: true,
                labelAlign: "center",
                color: "#70e5ea"
            }]
        },
        legend: {
            cursor: "pointer",
            itemclick: toggleDataSeries
        },
        data: [{
            type: "boxAndWhisker",
            toolTipContent: "<span style=\"color:#6D78AD\">{label}:</span> <br><b>Maximum:</b> {y[3]}<br><b>Q3:</b> {y[2]}<br><b>Median:</b> {y[4]}<br><b>Q1:</b> {y[1]}<br><b>Minimum:</b> {y[0]} <br><u>Total Counts:</u><br> <b>Quarter 4:</b> {quarter4}<br><b>Quarter 3:</b> {quarter3}<br><b>Quarter 2:</b> {quarter2}<br><b>Quarter 1:</b> {quarter1}",
            dataPoints: nucleotidePercentagesData
        },
        {
            type: "scatter",
            name: "Outlier Values",
            toolTipContent: "<span style=\"color:#C0504E\">{name}</span>: {y}",
            showInLegend: true,
            dataPoints: outlierdata,
            markerSize: 5
        }
    ]





//        nucleotidePercentagesData
    });
    nucleotideChart.render();
    enable_all_buttons();
    function toggleDataSeries(e) {
        if (typeof (e.dataSeries.visible) === "undefined" || e.dataSeries.visible) {
            e.dataSeries.visible = false;
        } else {
            e.dataSeries.visible = true;
        }
        e.chart.render();
    }
};
