Welcome to LoopDB
=================

LoopDB is a module for creating and storing DNA parts for Loop Assembly. It comes with a prebuilt
schema for implementing original Loop Assembly protocol, but is completely general and allows
alternative schema definitions ( see *manual-example* ).

For a quick introduction please read the **Getting started** section. Have a look at the **Tables**
section for a brief description of the underlying data structure. Check out example folders with
some code to get you started. You can clone just `examples` branch, without the sorce code by


::

    git clone -b examples git@github.com:HaseloffLab/LoopDB.git

Refer to **Advanced** for getting hold of **SQLAlchemy** functions.

**Any other questions?** Please file an issue, or contact me at m.delmans@gmail.com

Looking for LoopDesigner ?
==========================

It's on another `branch <https://github.com/HaseloffLab/LoopDB/tree/loopdesigner>`_.

Getting started
===============

Installing LoopDB
'''''''''''''''''

LoopDB is available through PyPi, so you can just

::

    pip install loopdb

Setting up a database
'''''''''''''''''''''

LoopDB is based on `SQLAlchemy <http://www.sqlalchemy.org>`__, which allows to use most of the
common database back-ends. We recommend using `PostgreSQL <https://www.postgresql.org>`__. After you
have installed PostgreSQL, or any other database server, create a new database. With PostgreSQL it
can be done with

::

    createdb loopdb

Initialising LoopDB
'''''''''''''''''''

Create a loopDB instance by

.. code:: python

    from loopDB import LoopDB

    loopDB = LoopDB( 'postgresql:///loopdb', clean = True )

where the first argument is the address of your database, and ``clean`` flag determines whether the
database should be emptied before initialisation. Default is ``clean=False``. 

Setting up a schema
'''''''''''''''''''

Now it's time to create a schema for the Loop Assembly. This involves setting up the restriction
enzyme, restriction enzyme site, base sequence and backbone tables, which will be used for defying
and storing parts. Read more about data structure in the **Tables** section. You have two options to
crate a schema, either using a json file or manually. Please refer to the *schema-example* and
*manual-example*.

Adding a part
'''''''''''''

You can define a new Level 0 part by using ``LoopDB.addPart()`` method, followed by the
``LoopDB.commit()`` method

.. code:: python

    myPart = loopDB.addPart(name = "My Part", seq = "ATG...", backbone = "L0-CDS")
    loopDB.commit()

where backbone is a name of a backbone you want to use for your part. If you are adding more than
one part, you can commit only onece after you issued several ``addPart()`` commands. Alternatively
you can load the part from a GenBank file, whcih will preserve and store all the sequence
annotations.

.. code:: python

    genBankPart = loopDB.addPart(name = "My GB Part",
                        gbFile = "part.gb", backbone = "L0-CDS")

Assembling a part
'''''''''''''''''

After you have defined several Level 0 parts you can use them to run an assembly using the same
``addPart()`` method. The only difference is that instead of specifying the sequence you should pass
a ``children`` argument containing names of the parts for assembly

.. code:: python

    myPromoter      = loopDB.addPart(name = "My Promoter",
                            seq = "GGA...", backbone = "L0-Prom5")
    myCDS           = loopDB.addPart(name = "My CDS",
                            seq = "ATG...", backbone = "L0-CDS")
    myTerminator    = loopDB.addPart(name = "My Part",
                            seq = "TCT...", backbone = "L0-Term3")
    myGene          = loopDB.addPart(name = "My Gene",
                            children = [myPromoter, myCDS, myTerminator], backbone = "Ly1")
    loopDB.commit()

One you have several Level 1 parts you can use the same method to assemble Level 2 parts and so on.
Furthermore, you can mix and match different levels, as long as their backbones have compatible
restriction sites. Having said that, make sure that the parts you pass to ``addPart()`` method are
in right order and are assembled in a compatible backbone, i.e. the 3' overhang of each successive
part matches 5' overhang of the previous part; and 3' overhang of the first part and 5' overhang of
the last part match the corresponding sequences of the template backbone.

Retrieving parts
''''''''''''''''

You can retrieve an existing part by passing part name to ``LoopDB.getPart()`` method. You can also
use, ``getBaseSeq()``, ``getBackbone()``, etc. to retrieve existing records for every LoopDB table.

.. code:: python

        myPart = loopDB.getPart("MyPart")

Further, you can access part's children or retrieve part sequence, by using the following properties

.. code:: python

        myPart.seq
        >> 'ATGGT...'
        myPart.fullSeq
        >> 'GTAGCAT ATG... GCTGAT'
        myPart.children
        >> [<tables.Part object at 0x10d5c8b10>, <tables.Part object at 0x10d5c8d50>]

The difference between ``seq`` and ``fullSeq`` is that the first one will return only the actual
sequence of the part, while the second one will return the complete sequence, including that of the
backbone. Additionally you can use ``record`` and ``fullRecord`` properties to get partial or
complete `Biopython <http://biopython.org>`__ ``SeqRecord`` that will include all the annotations.

Tables
======

LoopDB creates several tables behid the scene, which are used to store part elements.

RE Table
''''''''

RE Table stores Restriction enzyme definitions.

::

    RE Table
        name:       Name of restriction enzyme
        seq:        Recognition sequence of the enzyme

RES Table
'''''''''

RES Table stores pairs of restriction enzyme overhangs, that will be further used to define
receivers and adapters for Base sequences and Backbones.

::

    RES Table
        name:       Name of a restriction site
        site5:      Sequence of a 5' overhang
        site3:      Sequence of a 3' overhang

BaseSeq Table
'''''''''''''

BaseSeq Table stores definitions of Base sequences, which conceptually are meta-backbones composed
of a backbone sequence and receiver overhangs. In original Loop Assembly schema there are two Base
sequences: Ly (level odd) and Lx (Level even).

::

    BaseSeq Table
        name:       Name of a Base sequence
        seq:        Sequence
        receiver:   RES corresponding to the receiver overhangs

Backbone Table
''''''''''''''

Backbone table stores definitions of the backbones, which conceptually are variants of the Base
sequences, defined by unique adapter overhangs. In original Loop Assembly schema these are Lx1 - Lx4
and Ly1 - Ly4.

::

    Backbone Table
        name:       Name of a Backbone
        baseSeq:    Corresponding BaseSeq
        adapter:    RES corresponding to the adapter overhangs.
        *seq:       Backbone sequence
        *record:    Corresponding SeqRecord

Part Table
''''''''''

Here LoopDB stores all the parts.

::

    Part
        name:       Name of a Part
        backbone:   Corresponding backbone
        seq:        Original part sequence (Only for Level 0 parts)
        children:   List of references to subparts (For Level 1 and higher)
        *level:     Level of the part, defines as maximum level of the children + 1.
        *partSeq:   Paert sequence : a recursive sum of all part's children sequences (seq), including overhangs.
        *fullSeq:   Same as *partSeq but with backbone sequence included.
        *record:    SeqRecord with *partSeq as a sequnce pluss all annotations from the supplied gb files.
        *fullRecord: Same as *record but with backbone sequence and annotation.

*Note: Asterisks \* denote properties that are not stored in the database, but are genrated on the
fly.* # Advanced LoopDB is based on SQLAlchemy, which offers an advanced database querying system (
apart from many more other things ). You can get SQLAlchemy ``Session`` either by using
``LoopDB.session`` or creating a new session via ``LoopDB.Session()`` method.

.. code:: python

    loopDB = LoopDB(...)
    session = loopDB.Session()
    parts = session.query(Part).filter( ... ).all()
    session.close()

    # OR

    parts = loopDB.session.query(Part).filter( ... ).all()

For more information on querying see `SQLAlchemy
tutorial <https://docs.sqlalchemy.org/en/latest/orm/tutorial.html#querying>`__.

You can also have a look in ``tables.py`` file to see the definition of the SQLAlchemy tables and
their methods.

Any questions?
==============

Please feel free to file an issue or contact me at m.delmans@gmail.com
