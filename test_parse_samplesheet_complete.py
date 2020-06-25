import unittest

import results as rst
import parse_samplesheet_complete 

DATABASE = "captures.sqlite"
SAMPLESHEET = "SampleSheet.csv"


class TestSamplesheetParser(unittest.TestCase):
    def test_get_header(self):
        end_header = parse_samplesheet_complete.get_header(SAMPLESHEET)
        self.assertEqual(end_header, 21)
        
    def test_samplesheet_to_sample_genesis(self):
        sample_genesis = parse_samplesheet_complete.samplesheet_to_sample_genesis(SAMPLESHEET)
        self.assertEqual(sample_genesis, rst.SAMPLE_GENESIS)        
        
    def test_parse_samplesheet_for_pipeline(self):
        samplesheet_for_pipeline = parse_samplesheet_complete.parse_samplesheet_for_pipeline(SAMPLESHEET, DATABASE)
        self.assertEqual(samplesheet_for_pipeline, rst.PIPELINE_INFO)       
    
    def test_get_file_locations(self):
        samplesheet_for_pipeline = parse_samplesheet_complete.parse_samplesheet_for_pipeline(SAMPLESHEET, DATABASE)
        info_with_files = parse_samplesheet_complete.get_file_locations(samplesheet_for_pipeline, ".")
        self.assertEqual(info_with_files, rst.PIPELINE_INFO_FILES)

if __name__ == '__main__':
    unittest.main()
