var app = angular.module('app', ["angularplasmid"]);
// Globalizing seq
seq = "";
// Globalizing features
features = "";
socket = io();
downloadFrame = document.createElement('iframe');

downloadPart = function(){
	selected = w2ui["SideBar"].selected;
	if (selected != ''){
		downloadFrame.src = "/export?dbid=" + selected;
	}
}

$(function () {
	downloadFrame.style.display = 'none';
	document.body.appendChild(downloadFrame);

	$('#layout').w2layout({
		name: 'layout',
		panels: [
			{ type: 'left', size: 150, resizable: true, content: $().w2sidebar({
					name: "SideBar",
					nodes: [],
					parts: [],
					onClick: function(event){ renderPart(event.target) }	
				}) 
			},
			{ type: 'main' , toolbar: {
				items :[
					{ id: 'addNewButton', type: 'menu', caption: 'Add', img: 'w2ui-icon-plus', 
						items:[
							{"text" : "Add part", icon: "fa fa-circle-o-notch"},
							{"text" : "Assemble part", icon: "fa-infinity"}
						]
					}
				],
				onClick: function(event) {
					switch(event.target){
						case "addNewButton:Add part":
							renderAddNewForm();
							break;
						case "addNewButton:Assemble part":
							renderAssembleForm();
							break;
					}
				}
			}},
			{ type: 'preview', size: '50%', resizable: true}
		]
	});
	$("#layout_layout_panel_preview").attr("ng-controller", "preview");
	$("#layout_layout_panel_preview .w2ui-panel-content").attr("ng-include", "url");

	$("#layout_layout_panel_main").attr("ng-controller", "main");
	$("#layout_layout_panel_main .w2ui-panel-content").attr("ng-include", "url");

	console.log(w2ui);

	renderPartList();


	angular.element( function() {
		angular.bootstrap( document, ['app'] );
	});
});