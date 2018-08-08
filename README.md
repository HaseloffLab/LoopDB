# LoopDesigner

LoopDesigner is a web app for interacting with LoopDB. Currently, the following features are available:

* Adding a new part from a nucleotide sequence.
* Adding multiple parts from a Multi-FASTA file.
* Performing _in-silico_ loop assembly.
* Browsing plasmid maps of existing parts and constructs.
* Searching for a part by it's name.
* Automatic assembly protocol generation.

You can try a public version at [loopdesigner.herokuapp.com](http://loopdesigner.herokuapp.com)

# Installation

This section will explain how to install a local LoopDesigner server.

## Gettting Started

### Installing LoopDB

LoopDB is available through PyPi, so you can just

``` bash
pip install loopdb
```

### Setting up a database

LoopDB is based on [SQLAlchemy] (http://www.sqlalchemy.org), which allows to use most of the
common database back-ends. We recommend using [PostgreSQL](https://www.postgresql.org). After you
have installed PostgreSQL, or any other database server, create a new database. With PostgreSQL it
can be done with

``` bash
createdb loopdb
```

### Install LoopDesigner

Check out the `loopdesigner` branch:

``` bash
git clone https://github.com/HaseloffLab/LoopDB.git -b loopdesigner LoopDesigner
```

This will clone the files to LoopDesigner folder. Now install the requirements:

``` bash
cd LoopDesigner
pip install -r requirements.txt
```

### Starting the server

Now you should be ready to start the LoopDesigning server by running

```bash
python server.py
```

If you didn't get any errors you should be able to access the LoopDesigner by navigating to

http://127.0.0.1:8000

in your browser.

### Setting up a schema

The easiest way to define a schema is by using a JSON file. Have a look at the `schema.json` file for an example. The schema file should define four collections:

* RE: Restriction enzymes used for an assembly, e.g:

```javascript
{"name" : "BsaI", "seq" : "GGTCTCA"}
```

* RES: Pairs of restriction enzyme overhang sequences, which define either acceptors or recceivers, e.g.

```javascript
{"name": "AF", "re": "BsaI", "site5": "GGAG" , "site3": "CGCT"},
{"name": "ab", "re": "SapI", "site5": "ATG" , "site3": "GCA" }
```

* BaseSeq: Plasmid sequence with one of the receiver sequence pairs, defined in RES section. This corresponds to different assembly levels, e.g. for an odd level (pOdd):

```javascript
{"name": "pOdd-1", "receiver": "AF", "gbFile": "gb/pOdd-1.gb"}
```

where `gbFile` is the path to the plasmid sequence in the GenBank format.

* Backbone: A baseseq toogether with an adapter sequence, This corresponds to a particular vector, e.g. for pOdd-1

``` javascript
{"name": "pOdd-1", "baseSeq" : "pOdd", "adapter" : "ab" }
```

Once the schema file is complete, the schema can be applied to the database through an interactive python session:

```python
>> from loopDB import LoopDB # Importing the LoopDB module

>> loopDB = LoopDB( 'postgresql:///loopdb', clean = True) # Establishing the connection to the database
>> loopDB.initFromFile('schema.json') # Initialising from the schema file
```
