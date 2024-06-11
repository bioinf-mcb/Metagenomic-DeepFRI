from . import (test_alignment, test_bio_utils, test_install, test_mmseqs,
               test_pdb, test_predict, test_utils)


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_alignment))
    suite.addTests(loader.loadTestsFromModule(test_bio_utils))
    suite.addTests(loader.loadTestsFromModule(test_install))
    suite.addTests(loader.loadTestsFromModule(test_mmseqs))
    suite.addTests(loader.loadTestsFromModule(test_pdb))
    suite.addTests(loader.loadTestsFromModule(test_predict))
    suite.addTests(loader.loadTestsFromModule(test_utils))
    return suite
