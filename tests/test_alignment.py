import pytest

from meta_deepFRI.utils.search_alignments import alignment_sequences_identity


@pytest.mark.parametrize('query,target,expected', [('AASDS', 'AASDS', 1),
                                                   ('AASDS', 'ADSDS', 0.8),
                                                   ('AASSS', 'A-S-S', 0.6),
                                                   ('A--SS', 'A-S--', 0.2)])
def test_alignment_sequences_identity(query, target, expected):
    assert alignment_sequences_identity(query, target) == expected
