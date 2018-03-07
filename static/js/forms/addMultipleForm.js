function uploadFASTAForm(backbones){
	$().w2field('addType', 'label', function(options){});

	return $().w2form({
		name : "uploadFastaForm",
		header : "Add multiple parts",
		fields: [
			{ name: 'FASTA file', type: 'file', required: true,
				options: {
					max : 1
				},
				html:{
					attr: 'style=width:300px',
					column: 5,
					span: 3
				}
			},
		],
		actions: {
			"Upload" : function(event){
				fasta = atob(this.record["FASTA file"][0].content);
				var Fasta = require('biojs-io-fasta');
				seqs = Fasta.parse(fasta);
				
				if ('addMultipleForm' in w2ui){
					w2ui['addMultipleForm'].destroy();
				}

				w2ui['layout'].content('main', addMultipleForm(backbones, seqs) )
			}
		}
	})
}

addMultipleForm = function(backbones, seqs){
	fields = []	
	for (i=0; i<seqs.length; i++){
		this.fields.push({
			field: "label_" + i.toString(), idx: i, type: 'label', required: false,
			html: {
				column: 1,
				caption: seqs[i].name,
				attr: 'style=width:0px;visibility:hidden'
			}
		});

		this.fields.push({
			field: "baseSeq_" + i.toString(), idx: i, type: 'list', required: true,
			options: { 	items: backbones,
						renderDrop : function(item){
							return item.text
						}, 
			},
			html: {
				column: 2,
				caption: "BaseSeq"

			}
		});

		this.fields.push({
			field: "backbone_" + i.toString(), idx: i, type: 'list', required: true,
			options: { 	items: [],
						renderDrop : function(item){
							return item.text
						}, 
			},
			html: {
				column: 3,
				caption: "Backbone"
			}
		});
	}



	return $().w2form({
		name: "addMultipleForm",
		header: "Assign backbone",
		fields: fields,
		onChange: function(event){
			if (event.target.startsWith('baseSeq')){
				console.log( event );
				index = parseInt(event.target.split('_')[1]);
				this.fields[index*3 + 2].options.items = backbones.find(({text})=>(text == event.value_new.text)).backbones;
				this.refresh("backbone_"+index.toString());
			}
		},
		actions:{
			save : function(){


				if (this.validate(true).length == 0){
					for(i=0; i<this.fields.length / 3; i++){
						record={
							Backbone: this.record["backbone_"+i.toString()],
							Sequence: seqs[i].seq,
							"Part name": seqs[i].name
						};
						console.log("Record: ", record);
						socket.emit("addL0", record, function(response){
							if (response[0] == "OK"){
							}
							else{
							}
						});
					}
					renderPartList();
					renderPart(null);
				}
			}
		}
	});
}

