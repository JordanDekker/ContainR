function make_graphs(response) {
	var organism_name = response.node_name
	var dataLengthPoints = [];
	var dataIdentityPoints = [];
	var dataBitPoints = [];
	var Length_chart = new CanvasJS.Chart("Length", {
		animationEnabled: true,
		exportEnabled: true,
		zoomEnabled: true,
		theme: "light2",
		title: {
			text: "Distribution of sequence length (" + organism_name + ")"
		},
		axisX: {
			title: "Sequence length",
			interval: 10
		},
		axisY: {
			title: "Amount of sequences",
		},
		legend: {
			cursor: "pointer",
			dockInsidePlotArea: false
		},
		height: 300,
		width: 500,
		data: [{
			type: "column",
			indexLabel: "{y}",
			indexLabelFontColor: "#000000",
			indexLabelPlacement: "inside",
			dataPoints: dataLengthPoints
			}]
	});

	function add_length_Data(Length_data) {
		var g_length_data = Length_data.data;
		var arrayLength = g_length_data.length
		for (var i = 0; i < arrayLength; i++) {
			dataLengthPoints.push({
				x: g_length_data[i].length,
				y: g_length_data[i].values
			});
		}
		Length_chart.render();
	}

	var identity_chart = new CanvasJS.Chart("identity", {
		colorSet: "greenShades",
		animationEnabled: true,
		exportEnabled: true,
		theme: "light2",
		zoomEnabled: true,
		title: {
			text: "Distribution of identity (" + organism_name + ")"
		},
		axisY: {
			title: "Amount of sequences",
			prefix: "",
		},
		axisX: {
			title: "Percentage of similarities (%)",
			valueFormatString: "###",
			maximum: 100,
			minimum: 0
		},
		dataPointMaxWidth: 20,
		height: 300,
		width: 500,
		data: [{
			type: "column",
			indexLabel: "{y}",
			indexLabelOrientation: "vertical",
			indexLabelFontColor: "black",
			xValueFormatString: ("#'%'"),
			dataPoints: dataIdentityPoints
		}]
	});

	function add_identity_Data(identity_data) {
		var g_identity_data = identity_data.data;
		var arrayLength = g_identity_data.length
		for (var i = 0; i < arrayLength; i++) {
			dataIdentityPoints.push({
				x: g_identity_data[i].identity,
				y: g_identity_data[i].values
			});
		}
		identity_chart.render();
	}

	var bit_score_chart = new CanvasJS.Chart("Bit_score", {
		animationEnabled: true,
		zoomEnabled: true,
		exportEnabled: true,
		theme: "light2",
		title: {
			text: "Distribution of bit score (" + organism_name + ")"
		},
		axisY: {
			title: "Amount of sequences",
			includeZero: false,
			minimum: 0
		},
		axisX: {
			title: "Bit score",
		},
		height: 300,
		width: 500,
		data: [{
			type: "line",
			indexLabel: "{y}",
			dataPoints: dataBitPoints
			}]
	});

	function add_bit_Data(bit_data) {
		var g_bit_data = bit_data.data;
		var arrayLength = g_bit_data.length;
		for (var i = 0; i < arrayLength; i++) {
			dataBitPoints.push({
				x: g_bit_data[i].bitscore,
				y: g_bit_data[i].values
			});
		}
		bit_score_chart.render();
	}

	function add_taxonomy_links(tax_data) {
		// if tax id = none display message + link to tax db
		// if tax id = id: generate links to pages

		var ncbi_term = '<a href="https://www.ncbi.nlm.nih.gov/search/all/?term=' + organism_name.replace(" ", "%20") + '" target="_blank">NCBI</a>';
		var taxonomy_term = '<a href="https://www.ncbi.nlm.nih.gov/taxonomy/?term=' + organism_name.replace(" ", "+") + '" target="_blank">Taxonomy Database</a>';
		var pubmed_term = '<a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=' + organism_name.replace(" ", "+") + '" target="_blank">PubMed</a>';
		var ncbi_search = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/' + tax_data + '" target="_blank">NCBI search</a>';

		if (tax_data == "None" || tax_data == "nan") {
			document.getElementById("results-links").innerHTML = `
			<h2>Search links</h2>
			The selected node contains no reads with accession code. You can use the following links to search:<br>
			${ncbi_term}<br>
			${taxonomy_term}<br>
			${pubmed_term}
			`
		} else if (tax_data == "More") {
			document.getElementById("results-links").innerHTML = `
			<h2>Search links</h2>
			Not all reads in the selected node have the same accession code. You can use the following links to search:<br>
			${ncbi_term}<br>
			${taxonomy_term}<br>
			${pubmed_term}
			`
		} else if (tax_data == "unknown") {
			document.getElementById("results-links").innerHTML = `
			<h2>Search links</h2>
			All reads in the selected node have an unknown accession code.
			`
		} else {
			document.getElementById("results-links").innerHTML = `
			<h2>Search links</h2>
			Use the links below to search using the accession code:<br>
			${ncbi_search}<br>
			${taxonomy_term}<br>
			${pubmed_term}
			`
		}
	}

	function add_selection_options(response) {
		fileContent = "";
		for (i=0;i<response.fw_headers.length;i++){
			fileContent += ">" + response.fw_headers[i].substring(1) + "/1\n"
			fileContent += response.fw_seqs[i] + "\n"
		}
		for (i=0;i<response.rv_headers.length;i++){
			fileContent += ">" + response.rv_headers[i].substring(1) + "/2\n"
			fileContent += response.rv_seqs[i]+ "\n"
		}
		document.getElementById("selection-options").innerHTML = "<h2>Selection options</h2>Download the reads from the last selected node as a FASTA file. Forward reads have '/1' appended to the header, reverse reads '/2'.<br><a href=\"data:application/octet-stream;charset=utf-8;base64," + btoa(fileContent) + "\" download=\"export.fasta\"><button>Export selection</button></a>";
		document.getElementById("selection-options").innerHTML += " <button class=\"action_button\" onclick=\"delete_node()\">Remove selection</button><button id=\"remake_sunburst_button\" onclick=\"remake_sunburst_image()\">Remake sunburst</button>";
	}

	function split_json(response) {
		add_bit_Data(JSON.parse(response.count_bit));
		add_identity_Data(JSON.parse(response.count_id));
		add_length_Data(JSON.parse(response.count_length));
		add_selection_options(response);
		add_taxonomy_links(response.tax_id);
	}

	split_json(response);
};
