renderPart = function(partDBID){
	part = w2ui['SideBar'].parts.find( ({dbid}) => ( dbid == partDBID ) );
	angular.element( $("#layout_layout_panel_preview") ).scope().setPart( part );
	angular.element( $("#layout_layout_panel_preview") ).scope().$apply();
	// Added reset button to sequence-viewer
	// w2ui["layout"].content('main', "<div id='seqView'></div><button id='reset_btn' class='button right' onclick='reset();'>Reset</button>");

	w2ui["layout"].content('main', "<div id='seqView'></div>");

	socket.emit('getRecord', partDBID, function(record){
		// Assigning global seq
		seq = record["seq"];
		// Assigning global features
		features = record["features"];
		sequence = new Sequence(seq);
		sequence.render("#seqView", {"title" : part.name, "search" : true, "charsPerLine": 100, "sequenceMaxHeight": "200px"} );

		// Sorting features by position
		sort_features(features);
		// Highlighting features
		highlight(features);
		console.log(sequence);
		
	});
}

renderPartList = function(){
	socket.emit('getParts', '', function(parts){
		w2ui['SideBar'].parts = parts;
		w2ui['SideBar'].remove( w2ui['SideBar'].nodes.map( ({id}) => id ) );
		w2ui['SideBar'].add( parts.map( ({dbid, name}) => ({id: dbid, text: name}) )  );
	});
}

renderAddNewForm = function(){
	delete w2ui['addForm'];
	delete w2ui['addForm_tabs'];
	delete w2ui['addForm_toolbar'];
	socket.emit('getBackbones', function(backbones){
		console.log(backbones);
		var form = $().w2form({
			name : "addForm",
			header : "Add new L0 part",
			fields : [
				{ name: 'Part name', type: 'text', required: true,
					html:{
						attr: 'style=width:200px'
					}
				},
				{ name: 'GenBank file', type: 'file', required: true,
					options: {
						max : 1
					},
					html:{
						attr: 'style=width:200px',
					}
				},
				{ name: 'Backbone', type : 'list', required : true,
					options:{
						items : backbones,
						renderDrop : function(item){
							return item.text
						},
					},
					html:{
						attr: 'style=width:200px'
					}
				}
			],
			onValidate: function(event){
				if (w2ui['SideBar'].parts.find( ({name}) => name == this.record["Part name"] )){
					event.errors.push( {field: this.get("Part name"), error: "Part name already exists"} )
				}
			},
        	actions:{
        	
        		"Save" : function(event){
					console.log(this.record);
				
					if (this.validate(true).length == 0){
						socket.emit("addL0", this.record, function(response){
							if (response[0] == "OK"){
								w2ui['layout'].content('main', "<i class='fa fa-check-circle-o fa-5x'></i>");
								renderPartList();
							}
							else{
							}
						});
					}
					
	        	} 
        	}
		});
		w2ui['layout'].content('main', form);
		console.log(w2ui);
	});
}

assemblyForm = {
	init: function(){
		this.name  =  "assemblyForm";
		this.header = "Assemble new part";
		this.n =  0;
		this.part  =  {children : [], fullLength: 0, length: 0};
		this.tabs  = [
			{id : "tab0", caption : "Select Backbone"}
		];
		this.fields = [];
		this.page = 0;
		this.actions = {
			clear : {caption: "Clear", onClick : function(){
				assemblyForm.init();
			}},

			next : {caption : "Proceed", onClick : function(){

				socket.emit('getParts', this.record.selector.site3, function(parts){
					console.log(this.part);

					if (assemblyForm.n == 0 && parts.length == 0){
						this.onValidate = function(event){
							event.errors.push({ field: this.get("selector"), error : "Can't find parts to proceed assembly." });
						}
					}
					else{
						if (parts.length == 0 && this.record.selector.site3 != this.part.backbone.site5){
							this.onValidate = function(event){
								event.errors.push({ field: this.get("selector"), error: "Can't find parts to proceed assembly." });
							}
						}
						else{
							this.onValidate = function(event){};
						}
					}

					if (this.validate(true).length == 0){
						assemblyForm.tabs[ assemblyForm.tabs.length-1 ].caption = this.record.selector.text;
						if (assemblyForm.n == 0){
							assemblyForm.part.backbone = this.record.selector;
						}
						else{
							assemblyForm.part.children.push( this.record.selector );
							assemblyForm.part.length += this.record.selector.length;
						}
						
						assemblyForm.part.fullLength += this.record.selector.length;

						angular.element( $("#layout_layout_panel_preview") ).scope().setPart( assemblyForm.part );
						angular.element( $("#layout_layout_panel_preview") ).scope().$apply();


						assemblyForm.n++;
						for (i=0; i<assemblyForm.n; i++){
							assemblyForm.tabs[i].disabled = true;
						}
						assemblyForm.page = assemblyForm.n;

						if (parts.length != 0){
							assemblyForm.fields[0] = {
								field: "selector", type: "list", required: true,
								options:{
									items: parts
								},
								html:{
									page: assemblyForm.n,
									caption : "Select part " + assemblyForm.n
								}
							}
							assemblyForm.tabs.push( {id: "tab" + assemblyForm.n, caption: "Select part "+assemblyForm.n} )
						}
						else{
							assemblyForm.actions.next.caption = "Finalise";
							assemblyForm.actions.next.onClick = function(){

								this.onValidate = function(event){
									if (w2ui['SideBar'].parts.find( ({name}) => name == this.record["Part name"] )){
										event.errors.push( {field: this.get("Part name"), error: "Part name already exists"} )
									}
								};

								if (this.validate(true).length == 0){
									assemblyForm.part.name = this.record["Part name"];
									console.log( assemblyForm.part )
									socket.emit("submitAssembly", assemblyForm.part, function(response){
										if (response[0] == "OK"){
											console.log("RESPONE OK");
											w2ui['layout'].content('main', "<i class='fa fa-check-circle-o fa-5x'></i>");
											renderPartList();
										}
									});
								}
							}
							assemblyForm.fields = [
								{ field: "Part name", type: "text", required: true,
									html:{
										page: assemblyForm.n
									}
								}
							];
							assemblyForm.tabs.push( { id: "tabF", caption: "Set part name" } );
						}
						delete w2ui['assemblyForm'];
						delete w2ui['assemblyForm_tabs'];
						delete w2ui['assemblyForm_toolbar'];

						w2ui['layout'].content('main', $().w2form(assemblyForm) );

					}
					else{
						console.log("INVALID");
					}

				}.bind(this));

			}}
		}

		delete w2ui['assemblyForm'];
		delete w2ui['assemblyForm_tabs'];
		delete w2ui['assemblyForm_toolbar'];

		socket.emit('getBackbones', function(backbones){
			this.fields[0] = {
				field: "selector", type: "list", required: true,
				options:{
					items: backbones.filter( ({site3}) => !!site3 )
				},
				html:{
					caption : "Select backbone"
				}
			}
			w2ui['layout'].content('main', $().w2form(this));
			console.log(w2ui['assemblyForm']);
		}.bind(this));

	}
}

renderAssembleForm = function(){
	assemblyForm.init();
}
