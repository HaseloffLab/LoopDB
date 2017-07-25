function assemblyForm(backboneList){
	console.log("ASSEMBLY FORM INIT");
	form = {
		init: function(backboneList){
			this.name  =  "assemblyForm";
			this.header = "Assemble new part";
			this.n =  0;
			this.part  =  {children : [], fullLength: 0, length: 0};
			this.tabs  = [
				{id : "tab0", caption : "Select Backbone"}
			];
			this.fields = [];
			this.page = 0;

			this.proceedFunction = function(){
				console.log("Proceed THIS: ", this);
				socket.emit('getParts', this.record.selector.site3, function(parts){
						console.log(this.part);

						if (form.n == 0 && parts.length == 0){
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
							form.tabs[ form.tabs.length-1 ].caption = this.record.selector.text;
			
							if (form.n == 0){
								form.part.backbone = this.record.selector;
							}
							else{
								form.part.children.push( this.record.selector );
								form.part.length += this.record.selector.length;
							}
							
							form.part.fullLength += this.record.selector.length;
							
							plasmidScope = angular.element( $("#layout_layout_panel_preview") ).scope();
							console.log("Form.part: ", form.part);
							plasmidScope.setPart( form.part );
							plasmidScope.$apply();


							form.n++;
							for (i=0; i<form.n; i++){
								form.tabs[i].disabled = true;
							}
							form.page = form.n;

							if (parts.length != 0){
								form.fields[0] = {
									field: "selector", type: "list", required: true,
									options:{
										items: parts
									},
									html:{
										page: form.n,
										caption : "Select part " + form.n,
										attr: 'style=width:300px'
									}
								}
								form.tabs.push( {id: "tab" + form.n, caption: "Part "+form.n} )
							}
							else{
								form.actions.next = {};
								form.actions.next.caption = "Finalise";

								form.actions.next.onClick = this.finaliseFunction;
								form.fields = [
									{ field: "Part name", type: "text", required: true,
										html:{
											page: form.n,
											attr: 'style=width:200px'
										}
									}
								];
								form.tabs.push( { id: "tabF", caption: "Set part name" } );
							}

							w2ui['assemblyForm'].destroy();
							w2form = $().w2form(form);

							if (w2form.fields.find( ({field}) => field == "Part name" ) != null){
								console.log("Found FIELD");
								
								defName = w2form.part.children.map( ({text}) => text ).join("-");
								console.log(w2form.part.children);
								console.log(defName);
								w2form.record = {
									"Part name" : defName
								};
								w2form.refresh();
							}

							w2form.on({type: 'change', execute: 'after'}, w2form.proceedFunction);
							w2ui['layout'].content('main', w2form );
						}
						else{
							console.log("INVALID");
						}

					}.bind(this));
			},
			this.finaliseFunction = function(){
				this.onValidate = function(event){
					if (w2ui['sideBar'].parts.find( ({name}) => name == this.record["Part name"] )){
						event.errors.push( {field: this.get("Part name"), error: "Part name already exists"} )
					}
				};

				if (this.validate(true).length == 0){
					form.part.name = this.record["Part name"];
					console.log( form.part )
					socket.emit("submitAssembly", form.part, function(response){
						if (response[0] == "OK"){
							console.log("RESPONE OK");
							renderPartList(function(){
								w2ui['sideBar'].expandParents(response[1].dbid);
								w2ui['sideBar'].click(response[1].dbid);
							});
						}
					});
				}
			};

			this.actions = {
				back : {caption: "Back", onClick : function(){
					delete form.actions.next;

					if (form.n > 0){
						form.tabs.pop();
						
						form.n --;
						form.page = form.n;
						form.tabs[form.n].disabled = false;
						var site3 = null;

						if (form.n > 0){
							form.tabs[form.n].caption = "Part " + form.n;

							form.part.length -= form.part.children[form.part.children.length-1].length;
							form.part.fullLength -= form.part.children[form.part.children.length-1].length;
							form.part.children.pop();

							if (form.part.children.length > 0){
							site3 = form.part.children[form.part.children.length-1].site3;
							}
							else{
								site3 = form.part.backbone.site3;
							}

							socket.emit('getParts', site3, function(parts){
								form.fields[0] = {
										field: "selector", type: "list", required: true,
										options:{
											items: parts
										},
										html:{
											page: form.n,
											caption : "Select part " + form.n
										}
								}

								console.log("Back form: ", form);
								
								w2ui['assemblyForm'].destroy();
								w2form = $().w2form(form);
								w2form.on({type: 'change', execute: 'after'}, w2form.proceedFunction);
								
								w2ui['layout'].content('main', w2form );
								
								plasmidScope.setPart( form.part );
								plasmidScope.$apply();
							});
						}
						else{
							form.tabs[form.n].caption = "Select Backbone";
							form.fields[0] = {
								field: "selector", type: "list", required: true,
								options:{
									items: backbones.filter( ({site3}) => !!site3 )
								},
								html:{
									caption : "Backbone"
								}
							}
							w2ui['assemblyForm'].destroy();
							w2form = $().w2form(form);
							w2form.on({type: 'change', execute: 'after'}, w2form.proceedFunction);
							
							w2ui['layout'].content('main', w2form );
							
							plasmidScope.setPart( form.part );
							plasmidScope.$apply();
						}
					}
				}}
			} // actions
			
			console.log("W2UI: ", w2ui);

			if ('assemblyForm' in w2ui){
				w2ui['assemblyForm'].destroy();
			}

			console.log("Assembly form: ", form);
			backbones = [];
			for(i in backboneList){
				backbones = backbones.concat(backboneList[i].backbones);
			}

			form.fields[0] = {
				field: "selector", type: "list", required: true,
				options:{
					items: backbones.filter( ({site3}) => !!site3 )
				},
				html:{
					caption : "Backbone",
				}
			}
		} // init
	} // form

	form.init(backboneList);
	w2form = $().w2form(form);
	w2form.on({type: 'change', execute: 'after'}, w2form.proceedFunction);

	return w2form;
} // function