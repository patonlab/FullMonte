import unittest
import tempfile
import os
from FMMC import SDFWriter

class TestSDFWriter(unittest.TestCase):
    def test_basics(self):
        # Get a new temporary directory for each test in this class     
        tempdir      = tempfile.mkdtemp()
        # Assert that it's empty
        self.assertEqual(0, len(os.listdir(tempdir)))
        tempfilepath = os.path.join(tempdir, 'output.txt')

        self.assertFalse(os.path.isfile(tempfilepath))
        sdfwriter = SDFWriter(tempfilepath)
        self.assertTrue(os.path.isfile(tempfilepath))

        self.assertEqual(0, sdfwriter.num_individual_files)
        self.assertEqual(1, sdfwriter.counter)

        sdfwriter.write('hello')

        sdfwriter.next_conformation()
        self.assertEqual(2, sdfwriter.counter)

        # check 1 file in output dir
        self.assertEqual(1, len(os.listdir(tempdir)))
        sdfwriter.close()

    def test_individual_files(self):
        for num_individual_files in (1,2,3,4,5):
            # Get a fresh temp dir
            tempdir = tempfile.mkdtemp()
            # Assert that it's empty
            self.assertEqual(0, len(os.listdir(tempdir)))
            tempfilepath = os.path.join(tempdir, 'output.txt')
            self.assertFalse(os.path.isfile(tempfilepath))

            sdfwriter = SDFWriter(tempfilepath, num_individual_files=num_individual_files)
            self.assertTrue(os.path.isfile(tempfilepath))

            self.assertEqual(num_individual_files, sdfwriter.num_individual_files)

            for j in range(1,11):
                self.assertEqual(j, sdfwriter.counter)
                # check num files in output dir
                num_files = len(os.listdir(tempdir))
                if j > num_individual_files:
                    self.assertEqual(1+num_individual_files, num_files)
                else:
                    self.assertEqual(1+j, num_files)

                sdfwriter.write('hello%d' % j)
                sdfwriter.next_conformation()

            sdfwriter.close()

            # file 1 should contain hello1 and so on
            for j in range(1,num_individual_files):
                with open(os.path.join(tempdir, 'output_%d.txt' % j), 'r') as fd:
                    self.assertEqual('hello%d' % j, fd.read())
