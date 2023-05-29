import pytest

from mDeepFRI.alignment import insert_gaps
from mDeepFRI.alignment_utils import alignment_identity


@pytest.mark.parametrize('query,target,expected', [('AASDS', 'AASDS', 1),
                                                   ('AASDS', 'ADSDS', 0.8),
                                                   ('AASSS', 'A-S-S', 0.6),
                                                   ('A--SS', 'A-S--', 0.2)])
def test_alignment_identity(query, target, expected):
    seq_id = alignment_identity(query, target)
    assert round(seq_id, 6) == expected


@pytest.mark.parametrize(
    'query,reference,alignment,gapped_query,gapped_reference',
    [('AACT', 'AAT', 'MMDM', 'AACT', 'AA-T'),
     ('AAT', 'AATC', 'MMMI', 'AAT-', 'AATC'),
     ('AAT', 'FGTC', 'XXMI', 'AAT-', 'FGTC')])
def test_insert_gaps_full_length(query, reference, alignment, gapped_query,
                                 gapped_reference):
    gapped_seq, gapped_ref = insert_gaps(query, reference, alignment)
    assert gapped_seq == gapped_query
    assert gapped_ref == gapped_reference
    assert len(gapped_seq) == len(gapped_ref)
