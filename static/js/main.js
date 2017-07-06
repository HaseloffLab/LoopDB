// Initialising plasmidPreviewController from plasmid.js
angularPlasmid = angular.module("angularPlasmid", ["angularplasmid"]);
angularPlasmid.controller("plasmidPreview", plasmidPreviewController);

socket = io();

formConstructor = null;
downloadFrame = document.createElement('iframe');

downloadPart = function(){
	part = w2ui['sideBar'].parts.find( ({dbid}) => ( dbid == w2ui["sideBar"].selected ) );
	console.log(part);
	if ( part != undefined ){
		downloadFrame.src = "/export?dbid=" + part.dbid;
	}
}

deletePart = function(){
	part = w2ui['sideBar'].parts.find( ({dbid}) => ( dbid == w2ui["sideBar"].selected ) );
	console.log(part);
	w2confirm({
		msg: 'Are you sure you want to delete ' + part.text + '?',
		title: 'DNALooper'
	})
	.yes(function(){
		socket.emit('deletePart', part.dbid, function(response){
			renderPartList();
			renderPart(null);
		});
	})
	.no(function(){});
}

$(function () {

	// Adding download support
	downloadFrame.style.display = 'none';
	document.body.appendChild(downloadFrame);

	// Rendering layout defined in layout.js
	$('#layout').w2render('layout');
	
	// Setting bindings for angularPlasmid
	$( w2ui['layout'].el('preview') ).attr("ng-include", "url");

	angular.element( function() {
		angular.bootstrap( document, ['angularPlasmid'] );
	});

	// Getting Schema backbones
	socket.emit('getBackbones', function(backbones){
		formConstructor = new FormConstructor(backbones);
		console.log(formConstructor);
		renderPartList();
	});
})