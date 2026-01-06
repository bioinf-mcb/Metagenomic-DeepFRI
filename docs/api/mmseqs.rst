MMSeqs
======

.. currentmodule:: mDeepFRI.mmseqs

.. automodule:: mDeepFRI.mmseqs

.. autoclass:: mDeepFRI.mmseqs.QueryFile
   :special-members: __init__
   :members:

.. autoclass:: mDeepFRI.mmseqs.MMSeqsSearchResult
   :special-members: __init__
   :members:

Description
-----------

Python bindings and utilities for MMseqs2 fast similarity search. Includes
helpers for FoldComp databases, query filtering, and parsed search results.

Key Features
------------

- **Fast search**: Many-vs-many searches via MMseqs2 binaries bundled with the package
- **FoldComp support**: Extract FASTA sequences from compressed structure databases
- **Query handling**: Filter sequences by length and quality before search
- **Result parsing**: Access top hits, targets, and coverage/identity metrics

Example
-------

.. code-block:: python

   from mDeepFRI.mmseqs import QueryFile

   # Prepare query sequences
   qf = QueryFile("proteins.fasta")
   qf.load_sequences()
   qf.filter_sequences(min_length=30)

   # Search against an MMseqs2 database
   result = qf.search(
      database="pdb100.mmseqsDB",
      eval=1e-3,
      mmseqs_sensitivity=4,
      threads=8,
   )

   # Save raw results
   result.save("database_search/pdb100_results.tsv")
