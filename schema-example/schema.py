from loopDB import LoopDB

loopDB = LoopDB( 'postgresql:///loopdb', clean = True)
loopDB.initFromFile('schema.json')

cds = loopDB.addPart(name = "my CDS", backbone = "L0-CDS", seq = "ATGTAA")
loopDB.commit()

print "{0}: {1}".format( cds.dbid, cds.name )