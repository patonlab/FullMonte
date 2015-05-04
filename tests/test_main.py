from __future__ import print_function
import unittest
import tempfile
import os
import shutil
import FMMC
import logging

fmmc_logger = logging.getLogger('FMMC')
fmmc_logger.addHandler(logging.NullHandler())

class TestFMMCMain(unittest.TestCase):
    def setUp(self):
        self.hexane_path = os.path.join(os.path.dirname(__file__), 'fixtures', 'hexane.mol')
        self.tempdir = tempfile.mkdtemp()
        shutil.copy(self.hexane_path, self.tempdir)
        self.inputfile_path = os.path.join(self.tempdir, 'hexane.mol')

    def tearDown(self):
        pass

    def test_main(self):
        self.assertTrue(self.hexane_path.endswith('.mol'))
        self.assertEqual(['hexane.mol'], os.listdir(self.tempdir))

        # Run fmmc(). It should save results locally, into same directory as input file.
        # reproduce params as calculated in FMMC.py. NB imperfect. What if file contains >1 dot?
        filein   = self.inputfile_path.split(".")[0]
        filetype = self.inputfile_path.split(".")[1]
        FMMC.main(filein, filetype, maxstep=1)

        # hexane.dat, hexane.mol, hexane_fm.sdf
        self.assertEqual(3, len(os.listdir(self.tempdir)))