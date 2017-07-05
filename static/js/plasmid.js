function getAnnotation(part, annotation){
	var n = part.children.length;
	
	console.log( "Annotating part with " + n + " children" );

	if ( n>0 ){
		for (var i=0; i<n; i++){
			console.log(" Doing child " + i);
			annotation.pos += part.children[i].site5.length;
			annotation = getAnnotation(part.children[i], annotation);
			console.log("Done with child " + i + " out of " + n);
			console.log("Pos = " + annotation.pos);
		}
		annotation.pos += part.children[n-1].site3.length;
	}
	else{

		annotation.coverage.push({
			start	: annotation.pos,
			end  	: annotation.pos + part.length,
			bgcolor : part.color,
			color 	: "black"
		});

		annotation.pos += part.length;

		annotation.legend.push({
			name	: part.text,
			color	: part.color
		});
	}
	return annotation;
}

function highlight(part){

	// sequenceCoverage = [];
	// legend = [];

	annotation = getAnnotation(part, {coverage: [], legend: [], pos: 0});	

	// sequenceCoverage.push({
	// 	start:		start,
	// 	end:		end,
	// 	bgcolor:	color,
	// 	color:		"black",
	// 	underscore:	false   
	// });
	
	// legend.push(
	// 	{name: label, color: color, underscore: false}
	// );
		
	sequence.coverage(annotation.coverage);
	sequence.addLegend(annotation.legend);
	
}

app.controller('preview', function($scope){
	$scope.name = "";
	$scope.layers = [];
	$scope.markers = [];
	$scope.range = [];
	$scope.width = 25;
	$scope.radius = 100;
	$scope.length = 0;
	$scope.url = '/static/html/plasmid.html'
	$scope.click = function(event, marker){
		part = w2ui['SideBar'].parts.find( ({text}) => ( text == marker.labels[0].text ) );
		w2ui['SideBar'].click(part.dbid);
	}
	$scope.setPart = function(part){
		$scope.layers = [];
		$scope.range = [0];
		$scope.markers = [];

		$scope.layers.push( [part] );
		$scope.length = part.fullLength;
		$scope.name = part.name;

		partLen = part.length;

		nChildren = part.children.length;
		while (nChildren != 0){
			nChildren = 0;
			layer = [];
			start = 0;
			lastLayer = $scope.layers[ $scope.layers.length-1 ];
			for (i=0; i < lastLayer.length; i++){
				p = lastLayer[i];
				for (j=0; j < p.children.length; j++){
					pp = p.children[j];
					layer.push( pp );
					nChildren ++;
				}
			}
			if (nChildren > 0){
				$scope.layers.push(layer);
				$scope.range.push( $scope.layers.length-1 );
			}
		}

		console.log("Layers: ", $scope.layers);

		for(i=0; i < $scope.layers.length; i++){
			level = [];
			layer = $scope.layers[i];
			start = 0;
			for(j=0; j < layer.length; j++){
				part = layer[j];
				level.push( { start: start - 0.5*partLen, end: start + part.length - 0.5*partLen, text: part.text, color: part.color } );
				start += part.length;
			}
			$scope.markers.push(level);
			console.log(level);
		}

	};
})