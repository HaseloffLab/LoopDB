from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC

replace = {
			"GGTCTC" 	: ["GGaCTC","GGTCcC","GGTaTC"],
			"GAGACC" 	: ["GAGACt","GAGAtC" ,"GAGgCC"],
			"GCTCTTC"	: ["GCcCTTC","GCTCgTC","GCTgTTC"],
			"GAAGAGC"	: ["GAgGAGC","GAAGgGC","GAAaAGC"]
}

colorMap =['#99ccff','#F96381','#ffff80','#ff9966','#ccb3ff','#4db8ff','#ffcc66','#d5ff80','#C6A49A','#56FAAF']

partColors = {		
	'prom5' 	: "#ffff9b",
	'prom'  	: "#ffffc8",
	'utr5'		: "#c8f0fa",
	
	'cds'		: "#a0ffa0",
	'cds-' 		: "#beffbe",
	'cds1'		: "#beffbe",
	'cds2'		: "#beffbe",
	'ctag'		: "#beffbe",
	
	'term'  	: "#fac8c8",
	'utr3'		: "#c8f0fa",
	'term3' 	: "#ffa5a5",


	'TU'		: "#beffbe",
	
	'ab'		: "#e54e51",
	'bg'		: "#f18a5e",
	'ge'		: "#fbc250",	
	'ew'		: "#fbf083",
	'AB'		: "#3765a5",
	'BC'		: "#8ca9cf",
	'CE'		: "#cce0f4",	
	'EF'		: "#acd8c8",
}

def recfind(pattern, string, start=0):
    pos = string.find(pattern, start)
    if pos == -1:
        return []

    return [pos] + recfind(pattern, string, pos + len(pattern))

def check(seq, sites):
	found = 0
	for site in sites:
		a = seq.find(site, 0, len(seq))
		if a != -1:
			found = found + 1
	return found

def domesticate(rec):
	sites = {}

	# Finding restriction sites
	for pattern in replace:
		sites[pattern] = recfind(pattern, rec.seq.upper())

	# Replacing and labeling restriction sites
	n = 0
	for pattern in sites:
		for site in sites[pattern]:
			rec = rec[:site] + Seq( replace[pattern][ (site-1) % 3 ], IUPAC.unambiguous_dna ) + rec[site + len(pattern):]
			
			i = n % 10
			n = n + 1

			recFeature = SeqFeature( FeatureLocation(site, site + len(pattern)),
				type = "misc_feature", id = "Restriction site",
					qualifiers = {"label" : ["Restriction site removed"], "ApEinfo_fwdcolor": [colorMap[i]] } )
			rec.features.append(recFeature)

	# Checking there are no more sites present
	if check(rec.seq.upper(), sites) > 0:
		return None

	return rec