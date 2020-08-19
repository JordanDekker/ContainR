//window.onload = function () {
function forwardReverseCompare(data){
    var dataPoints = [];
    CanvasJS.addColorSet("greenShades",             //creates a colorSet which colors the columns
                    [//colorSet Array

                    "#2F4F4F",
                    "#008080",
                    "#2E8B57",
                    "#3CB371",
                    "#90EE90"
                    ]);
    var pairsChart = new CanvasJS.Chart("identityImage", {
        colorSet: "greenShades",                                                //activates colorSet
        animationEnabled: true,                                                 //enables startup animation
        exportEnabled: true,
        theme: "light2", // "light1", "light2", "dark1", "dark2"                //sets the theme of the graph
        zoomEnabled: true,                                                      //enables zoom function
        title: {
            text: "Similarities of the Forward and Reverse complement" //title
        },
        axisY: {
            title: "Amount of sequences",           //title
        },
        axisX: {
            title: "Percentage of similarities (%)",    //title
            valueFormatString:"###",                //sets the x-axis values to only consist of 3 numbers max
            maximum: 100,                           //sets a limit on the maximum value of the x-axis
            minimum: 0                             //sets a limit on the minimum value of the y-axis

        },
        dataPointMaxWidth: 20,
        data: [{
            type: "column",                         //sets the kind of graph you will get
            indexLabel: "{y}",                      //gives the columns a label based on their Y-axis value
            indexLabelOrientation: "vertical",      //makes the above named label appear vertical
            indexLabelFontColor: "black",           //makes the above named label black
            xValueFormatString: ("#'%'"),           //when you go over a column with your mouse a % will show before the X-axis value
            dataPoints: dataPoints
        }]
    });
    /*
    Function takes json array (example: [
      {
        "percentage": 90,
        "amount": 100
      }]
      )
      it loops through the array and takes the *percentage* for the x-axis and the *amount* for the y-axis
    */
    function addData(data) {
        var g_data = data.data;
        var arrayLength = g_data.length
        for (var i = 0; i < arrayLength; i++) {
            dataPoints.push({
                x: g_data[i].perc,
                y: g_data[i].overlap_identity_perc
            });
        }
        pairsChart.render();

    }
    addData(data);
    enable_all_buttons();
}