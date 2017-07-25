// Load parts from server and render into SideBar

function sortedDict(dict){
	var keys = Object.keys(dict);
	keys.sort( (a,b) => a.localeCompare(b) );
	
	var sortedDict = {};

	for(var i=0; i < keys.length; i++){
		sortedDict[keys[i]] = dict[ keys[i] ];
	}

	return sortedDict;
}

function sortedParts(parts){
	
	var sortedParts = [];
	var sortedKeys = parts.map( ({text}) => text);
	sortedKeys.sort( (a,b) => a.localeCompare(b) );
	for(var i=0; i<sortedKeys.length; i++){
		key = sortedKeys[i];
		sortedParts.push( parts.find( ({text}) => text==key) );
	}

	return sortedParts;
}

renderPartList = function(callback){

	console.log("Rendering part list with callback", callback);
	socket.emit('getParts', '', function(parts){
		w2ui['sideBar'].parts = parts;
		w2ui['sideBar'].remove.apply(w2ui['sideBar'], w2ui['sideBar'].nodes.map( ({id}) => id ));
		
		$('input[type=combo]').w2field('combo', {
			items: parts.map( ({text, id}) => ({text:text, id: id}) )
		});
		$('input[type=combo]').on('change', function(event){
			part = w2ui['sideBar'].parts.find( ({text}) => ( text == event.currentTarget.value ) );
			if (part.dbid != w2ui['sideBar'].selected){
				w2ui['sideBar'].expandParents(part.dbid);
				w2ui['sideBar'].click(part.dbid);
				$('#searchField input')['0'].value="";
			}
		});

		levels = [];
		
		for(i=0; i<parts.length; i++){
			part = parts[i];
			level = part.level;
			
			if (levels[level] == undefined){
				levels[level] = {};
			}
			
			bName = part.backbone.name;
			
			if (!(bName in levels[level])){
				levels[level][bName] = [];
			}

			levels[level][bName].push({id: part.dbid, text: part.name, nodes:[], part:true, icon:"fa fa-circle-o-notch"});
		}

		for(i=0; i<levels.length; i++){
			if (levels[i] == undefined){
				levels[i] = {};
			}

			levelID = "Level " + i.toString()

			w2ui['sideBar'].add({id: levelID, text: levelID, group: true, part:false }); 

			for(var bName in sortedDict(levels[i]) ){
				bID = levelID + "-" + bName;
				w2ui['sideBar'].insert(levelID, null, {id: bID, text: bName, part:false, img: 'icon-folder'});
				w2ui['sideBar'].insert(bID, null, sortedParts(levels[i][bName]) );
			}
		}
		console.log("Doing Callback", callback);
		if (callback != null){
			callback();
		}
		
		console.log("Done Callback");
	});
}

function getAnnotation(part, annotation){
	var n = part.children.length;
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
	annotation = getAnnotation(part, {coverage: [], legend: [], pos: 0});	

	sequence.coverage(annotation.coverage);
	sequence.addLegend(annotation.legend);
}

// Rendering Plamsid map, and Sequence viewer
renderPart = function(part){
	plasmidScope = angular.element( $("#layout_layout_panel_preview") ).scope();
	console.log("Part: ", part);

	if (part == null){
		plasmidScope.clear();
		plasmidScope.$apply();
		w2ui["layout"].content("main", "");
	}

	else{
		// Rendering plasmid
		plasmidScope.setPart( part );
		plasmidScope.$apply();

		// Rendering Sequence Features and Table
		w2ui["layout"].content('main', "<div id='seqView'></div><div id='partTable'><table><tbody class='list'></tbody></table></div>");

		socket.emit('getRecord', part.dbid, function(record){
			seq = record.seq.toUpperCase();

			features = record.features;
			sequence = new Sequence(seq);

			sequence.render("#seqView", {"title" : part.name, "search" : true, "charsPerLine": 100, "sequenceMaxHeight": "400px"} );
			highlight(part);

			nChildren = part.children.length;

			if (nChildren != 0){
				valueNames = ["name", "backbone"];
				values = [{ name: part.name, backbone: part.backbone.name }, { name: "Length", backbone: part.backbone.length.toString() + " bp" },
							{name: "Concentration", backbone: Math.round(part.backbone.length / 200).toString() + " ng/uL" }];
				
				item = "<tr><td class='name'></td><td class='backbone'></td>"
				
				for(i=0; i<nChildren; i++){
					name = "Part " + (i+1).toString();
					valueNames.push(name);
					item += "<td class='" + name + "'></td>"
					values[0][name] = part["children"][i].name;
					values[1][name] = part["children"][i].fullLength.toString() + " bp";
					values[2][name] = Math.round(part["children"][i].fullLength / 100).toString() + " ng/uL";
				}

				item += "</tr>";

				options = { valueNames: valueNames, item: item };

				console.log(options, values);

				partList = new List("partTable", options, values);
			}
		});
	}
}