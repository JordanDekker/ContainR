// main sunburst function
function build_sunburst(json) {
    var width = 600;
    var height = 600;
    var radius = Math.min(width, height) / 2;

    // Total size of all segments; we set this later, after loading the data.
    var totalSize = 0;

    var vis = d3.select("#chart").append("svg:svg")
        .attr("width", width)
        .attr("height", height)
        .append("svg:g")
        .attr("id", "container")
        .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

    var partition = d3.partition()
        .size([2 * Math.PI, radius * radius]);

    var arc = d3.arc()
        .startAngle(function (d) { return d.x0; })
        .endAngle(function (d) { return d.x1; })
        .innerRadius(function (d) { return Math.sqrt(d.y0); })
        .outerRadius(function (d) { return Math.sqrt(d.y1); });

    // Use d3.text and d3.csvParseRows so that we do not need to have a header
    // row, and can receive the csv as an array of arrays.
    function load(json) {
        createVisualization(json);
    };

    // Main function to draw and set up the visualization, once we have the data.
    function createVisualization(json) {
        // Bounding circle underneath the sunburst, to make it easier to detect
        // when the mouse leaves the parent g.
        vis.append("svg:circle")
            .attr("r", radius)
            .style("opacity", 0);

        // Turn the data into a d3 hierarchy and calculate the sums.
        var root = d3.hierarchy(json)
            .sum(function (d) { return d.size; })
            .sort(function (a, b) { return b.value - a.value; });

        // For efficiency, filter nodes to keep only those large enough to see.
        var nodes = partition(root).descendants()
            .filter(function (d) {
                return (d.x1 - d.x0 > 0.005); // 0.005 radians = 0.29 degrees
            });

        var path = vis.data([json]).selectAll("path")
            .data(nodes)
            .enter().append("svg:path")
            .attr("display", function (d) { return d.depth ? null : "none"; })
            .attr("d", arc)
            .attr("fill-rule", "evenodd")
            .style("fill", function (d) { return "#" + CryptoJS.enc.Hex.stringify(CryptoJS.SHA1(d.data.name)).substring(0, 6); })
            .style("opacity", 1)
            .on("mouseover", mouseover);

        // Add the mouseleave handler to the bounding circle.
        d3.select("#container").on("mouseleave", mouseleave);

        // Get total size of the tree = value of root node from partition.
        totalSize = path.datum().value;
    };

    // Fade all but the current sequence, and show it in the breadcrumb trail.
    function mouseover(d) {
        var percentage = (100 * d.value / totalSize).toPrecision(3);
        var percentageString = percentage + "%";
        if (percentage < 0.1) {
            percentageString = "< 0.1%";
        }

        d3.select("#percentage")
            .text(percentageString);

        d3.select("#explanation")
            .style("visibility", "");

        var sequenceArray = d.ancestors().reverse();
        sequenceArray.shift(); // remove root node from the array
        updateBreadcrumbs(sequenceArray, percentageString);

        // Fade all the segments.
        d3.selectAll("path")
            .style("opacity", 0.3);

        // Then highlight only those that are an ancestor of the current segment.
        vis.selectAll("path")
            .filter(function (node) {
                return (sequenceArray.indexOf(node) >= 0);
            })
            .style("opacity", 1);

        // on click send a json string with the taxonomic data of the selected bacterium
        vis.on("click", function (d) {
            data_ir = sequenceArray
            var i
            var data_json = []
            for (i = 0; i < data_ir.length; i++) {
                data_json.push(data_ir[i]['data']["name"])
            }
            make_bacterial_vizualisation(data_json)

        });
    }

    // Restore everything to full opacity when moving off the visualization.
    function mouseleave(d) {
        // Deactivate all segments during transition.
        d3.selectAll("path").on("mouseover", null);

        // Transition each segment to full opacity and then reactivate it.
        d3.selectAll("path")
            .transition()
            .style("opacity", 1)
            .on("end", function () {
                d3.select(this).on("mouseover", mouseover);
            });

        d3.select("#explanation")
            .style("visibility", "hidden");
    }

    // Update the breadcrumb trail to show the current sequence and percentage.
    function updateBreadcrumbs(nodeArray, percentageString) {
        $("#breadcrumb").empty();

        nodeArray.forEach(function(item) {
            color = CryptoJS.enc.Hex.stringify(CryptoJS.SHA1(item.data.name)).substring(0, 6)
            contrast = Math.round(((parseInt(color.substring(0, 2), 16) * 299) + (parseInt(color.substring(2, 4), 16) * 587) + (parseInt(color.substring(4, 6), 16) * 114)) / 1000);

            $("#breadcrumb").append("<li style=\"color:#" + (contrast > 128 ? "000" : "fff") + ";background-color:#" + color + "\">" + item.data.name + "</li>");
        });
    }

    load(json);
};
