from loopDB import LoopDB

loopDB = LoopDB( 'postgresql:///loopdb', clean = True)
loopDB.initFromFile('schema.json')

prom5_35s = loopDB.addPart(name = '35s Promoter', 	gbFile = 'gb/35S_PROM5.ape', backbone = 'L0-Prom5')
cds_eGFP  = loopDB.addPart(name = 'eGFP',			gbFile = 'gb/eGFP_CDS.ape',  backbone = 'L0-CDS-')
tag_lti	  = loopDB.addPart(name = 'LTI Tag',		gbFile = 'gb/Lti6b-TAG.ape', backbone = 'L0-CTAG')
term_35s  = loopDB.addPart(name = '35s Terminator', gbFile = 'gb/35S_TERM.ape',  backbone = 'L0-TERM3')

loopDB.commit()
