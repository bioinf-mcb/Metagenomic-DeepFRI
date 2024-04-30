Examples
========

MMSeqs2 search
--------------

.. code-block:: python

    from mDeepFRI.mmseqs import QueryFile

    # Create a query file
    query = QueryFile("sequences.fasta")
    # Load sequences to manipulate on them
    query.load_sequences()
    # filter sequences under 30 amino acids
    query.filter_sequences(min_length = 30)
    # search against MMSeqs2 database
    result = query.search("mmseqs.db", eval=10e-3, sensitivity=4, threads=8)
    # save results
    result.save("output.tsv")
