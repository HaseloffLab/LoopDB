/*
 * msa-seqtools
 * https://github.com/greenify/msa-seqtools
 *
 * Copyright (c) 2014 greenify
 * Licensed under the MIT license.
 */

var st = {};
define(function (require, exports, module) {
	module.exports = st;
});

/****
 * Seems to be lots of different ways to format FASTA headers. 
 * 
 * Generally there's an ID and a DESCRIPTION
 *   >ID DESCRIPTION
 * 
 *   >(parts|of|ID) (DESCRIPTION with optional key=values)
 *   
 * This is complicated by the fact that the "values" in the description can have spaces
 * e.g. OS=Arabidopsis thaliana GN=CCD8
 * 
 ****
*/

// extract IDs and push them to the meta dict
st.getMeta = function(label) {
  
	var full_id = false, full_desc = false;
	var name, ids = {}, details = {}, description;

// 	console.log( "getMeta.label: ", label );
  
	var label_parts = label.split(" ");

	if ( label_parts.length >= 1 ) {
		full_id   = label_parts.shift();     // everything up to the first white space
		full_desc = label_parts.join(" ");   // everything else
	}
	else {
		full_id = label; 
	}
	
// 	console.log( "full_id", full_id );
// 	console.log( "full_desc", full_desc );
	
	if ( full_id ) {
		var id_parts = full_id.split('|');
		
		// the last item is the accession
		name = id_parts.pop(); 
		
		details.en = name;
		
		// everything else should be pairs: db|id
		while ( id_parts.length != 0 ) {
			var db = id_parts.shift();
			var id = id_parts.shift();
			ids[ db ] = id;
		}
	}
	else {
		name = full_id;
	}

	if ( full_desc ) {
	
		var kv_parts = full_desc.split('=');
	
		if ( kv_parts.length > 1 ) {
			
			var current_key, next_key;
			var kv;
			var kv_idx_max = kv_parts.length - 1;
			var kv_idx = 0;
			kv_parts.forEach( function( value_and_maybe_next_key ) {
				
				value_and_maybe_next_key = value_and_maybe_next_key.trim();
				
				var value_parts = value_and_maybe_next_key.split(" ");
				var value;
				if ( value_parts.length > 1 ) {
					next_key = value_parts.pop();
					value = value_parts.join(' ');
				}
				else {
					value = value_and_maybe_next_key;
				}

				if ( current_key ) {
					var key = current_key.toLowerCase();
					details[ key ] = value;
					//console.log( "details[" + key + "] = " + value );
				}
				else {
					description = value;
					//console.log( "description=" + value );
				}
				current_key = next_key;
			});
		}
		else {
			description = kv_parts.shift();
		}
	}
	
	var meta = {
		name: name,
		ids: ids,
		details: details,
	};
	
	if ( description ) {
		meta.desc = description
	}
	
// 	console.log( "meta", meta );
	
	return meta;
};

var findSepInArr = function(arr, sep) {
  for (var i = 0; i < arr.lenght; i++) {
    if (arr[i].indexOf(i)) {
      return i;
    }
  }
  return arr.length - 1;
};

var strToDict = function(str, sep, toJoin) {
  toJoin = toJoin || {};
  var entries = str.split(sep);
  toJoin[entries[0].toLowerCase()] = entries[1];
  return toJoin;
};

var identDB = {
  "sp": {
    link: "http://www.uniprot.org/%s",
    name: "Uniprot"
  },
  "tr": {
    link: "http://www.uniprot.org/%s",
    name: "Trembl"
  },
  "gb": {
    link: "http://www.ncbi.nlm.nih.gov/nuccore/%s",
    name: "Genbank"
  },
  "pdb": {
    link: "http://www.rcsb.org/pdb/explore/explore.do?structureId=%s",
    name: "PDB"
  }
};

st.buildLinks = function(meta) {
  var links = {};
  meta = meta || {};
  Object.keys(meta).forEach(function(id) {
    if (id in identDB) {
      var entry = identDB[id];
      var link = entry.link.replace("%s", meta[id]);
      links[entry.name] = link;
    }
  });
  return links;
};


// search for a text
st.contains = function(text, search) {
  return ''.indexOf.call(text, search, 0) !== -1;
};

// split after e.g. 80 chars
st.splitNChars = function(txt, num) {
  var i, _ref;
  num = num || 80;
  var result = [];
  for (i = 0, _ref = txt.length - 1; i <= _ref; i += num) {
    result.push(txt.substr(i, num));
  }
  return result;
};

st.reverse = function(seq) {
  return seq.split('').reverse().join('');
}

st.complement = function(seq) {
  var newSeq = seq + "";
  var replacements = [
    // cg
    [/g/g, "0"],
    [/c/g, "1"],
    [/0/g, "c"],
    [/1/g, "g"],
    // CG
    [/G/g, "0"],
    [/C/g, "1"],
    [/0/g, "C"],
    [/1/g, "G"],
    // at
    [/a/g, "0"],
    [/t/g, "1"],
    [/0/g, "t"],
    [/1/g, "a"],
    // AT
    [/A/g, "0"],
    [/T/g, "1"],
    [/0/g, "T"],
    [/1/g, "A"],
  ];

  for(var rep in replacements){
    newSeq = newSeq.replace(replacements[rep][0], replacements[rep][1]);
  }
  return newSeq;
}

st.reverseComplement = function(seq){
  return st.reverse(st.complement(seq));
}

st.model = function Seq(seq, name, id) {
  this.seq = seq;
  this.name = name;
  this.id = id;
  this.ids = {};
};