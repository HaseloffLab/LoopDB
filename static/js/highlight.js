function sort_features(features){
	features.sort(function(a, b) {
    	return a.start - b.start;
	});
}

sequenceCoverage = [];
legend = [];

function highlight(features){
	sort_features(features);
	sequenceCoverage = [];
	legend = [];
	var previous = 0;
	
	var arrayLength = features.length

	for (var i = 0; i < arrayLength; i++) {

		var start = features[i].start;
		var end = features[i].end;
		var label = features[i].label;

		if (start < previous){
			console.log("Skipping label "+label+" due to being inside another feature");
			continue;
		}

		previous = end;

		color = features[i].color;

		sequenceCoverage.push({
			start:		start,
			end:		end,
			bgcolor:	color,
			color:		"black",
			underscore:	false   
		});
		
		legend.push(
			{name: label, color: color, underscore: false}
		);
		
		$('.coverageLegend').empty();

		sequence.coverage(sequenceCoverage);
		sequence.addLegend(legend);
	}
}
// Do not know to make distinction from first highlighting to the selected highlighting for reset button to work
function highlight_selected(features){
	var sequenceCoverage = [];
	var legend = [];
	var subSeq = "";
	var previous = 0;
	
	var arrayLength = features.length

	for (var i = 0; i < arrayLength; i++) {

		var start = features[i].start;
		var end = features[i].end;
		var label = features[i].label;
		
		if (start < previous){
			console.log("Skipping label "+label+" due to being inside another feature");
			continue;
		}

		var index = sequenceCoverage.length % 10;

		previous = end;

		color = typeColorsSelected[index];

		sequenceCoverage.push({
			start:		start,
			end:		end,
			bgcolor:	color,
			color:		"black",
			underscore:	false
		});
		
		legend.push(
			{name: label, color: color, underscore: false}
		);

		subSeq = subSeq + seq.substring(start, end);
	
		sequence.coverage(sequenceCoverage)
		sequence.addLegend(legend);


	}
}

function reset(){
	// Dirty hack to remove old Legend. Problems when original features are == 0. I think when you assemble, no features are defined nor retained when storing in database.
	$('.coverageLegend').empty();
	sequence.coverage(sequenceCoverage);
	sequence.addLegend(legend);
	console.log(sequenceCoverage);
};	