"""Test annotation parsing."""
import os
import unittest

import pandas as pd
from imaps.base.constants.sheet import METHOD
from imaps.base.sheet import DescriptorSheet


class TestDescriptorSheet(unittest.TestCase):
    def setUp(self):
        tests_dir = os.path.dirname(os.path.realpath(__file__))
        self.example_annotation_path = os.path.join(tests_dir, "data/sample_annotation.xlsx")
        self.example_annotation_tab_path = os.path.join(tests_dir, "data/sample_annotation.tab.gz")
        self.example_annotation = DescriptorSheet(self.example_annotation_path)
        self.correct_content = [
            [
                "Sample_1",
                "Test collection Joe",
                "Blah blah",
                "Super Smarty",
                "The Boss",
                "iCLIP",
                "It was magic",
                "TARDBP",
                "HEK293",
                "",
                "Homo sapiens",
                "NNNN,GTAAC_0,NNNNN",
                "AGATCGGAAG_1,AGCGGTTCAG_2",
                "HiSeq",
                "mouse anti-Pseudouridine",
                "irCLIP_ddRT_42",
                "L3-GTC",
                "",
                "",
                "",
                "1995.0",
                "",
                "1.0",
                "",
                "122.0",
                "10.0",
                "2020-04-10",
                "",
            ],
            [
                "Sample_2",
                "Test collection Joe",
                "Blah blah",
                "Super Smarty",
                "The Boss",
                "iCLIP",
                "It was magic",
                "TARDBP-GFP",
                "CEM_SS",
                "",
                "Homo sapiens",
                "NNNN,GTAAC_0,NNNNN",
                "AGATCGGAAG_1,AGCGGTTCAG_2",
                "HiSeq",
                "mouse anti-Pseudouridine",
                "irCLIP_ddRT_72",
                "L3-GGA",
                "",
                "",
                "",
                "2000.0",
                "",
                "2.0",
                "",
                "123.0",
                "15.0",
                "2020-04-10",
                "",
            ],
            [
                "Sample_3",
                "Test collection Joe",
                "Blah blah",
                "Super Smarty",
                "The Boss",
                "iCLIP",
                "It was magic",
                "ctrl-pseudouridine",
                "Cal51",
                "Serious",
                "Homo sapiens",
                "NNNN,CCGGA_0,NNN",
                "AGATCGGAAG_1,AGCGGTTCAG_2",
                "HiSeq",
                "mouse anti-Pseudouridine",
                "irCLIP_ddRT_36",
                "L3",
                "yes",
                "",
                "",
                "2010.0",
                "",
                "3.0",
                "",
                "124.0",
                "20.0",
                "2020-04-10",
                "",
            ],
        ]
        self.df_correct_content = pd.DataFrame(self.correct_content, columns=self.example_annotation.column_names)

    def change_column_content(self, column_name, replacement):
        """Replace the content in a column."""
        if column_name not in self.example_annotation.column_names:
            print("{} is not present in column names.".format(column_name))

        df = pd.DataFrame(self.example_annotation.content, columns=self.example_annotation.column_names)
        df[column_name] = replacement
        self.example_annotation.content = df.values.tolist()

    def test_fname(self):
        """Test if the saved filename is correct."""
        self.assertEqual(self.example_annotation.fname, self.example_annotation_path)

    def test_extension(self):
        """Test if the saved format is correct."""
        self.assertEqual(self.example_annotation.extension, "xlsx")

        # Test if .gz gets recognized.
        self.example_annotation.fname = self.example_annotation_tab_path
        self.example_annotation.get_extension()

        self.assertEqual(self.example_annotation.extension, "tab.gz")

    def test_content(self):
        """Test if the content is correctly read."""
        self.assertEqual(self.example_annotation.content, self.correct_content)

        # Test if the content is correctly read from tab.gz.
        example_annotation_tab = DescriptorSheet(self.example_annotation_tab_path)

        self.assertEqual(example_annotation_tab.content, self.correct_content)
        self.assertEqual(example_annotation_tab.column_names, self.example_annotation.column_names)

    def test_get_column(self):
        """Test if the right column is returned."""
        col = self.example_annotation.get_column(column_name="Protein")
        self.assertEqual(col, self.df_correct_content["Protein"].tolist())

    def test_get_element(self):
        """Test if the right element is returned."""
        element = self.example_annotation.get_element(sample_name="Sample_1", column_name="Protein")
        self.assertEqual(element, "TARDBP")

    def test_get_barcodes3(self):
        """Test if the 3' barcodes and number of tolerated mismatches are returned."""
        self.assertEqual(self.example_annotation.get_barcodes3(), (["GTC", "GGA", ""], 3))

    def test_get_barcodes5(self):
        """Test if the 5' barcodes and number of tolerated mismatches are returned."""
        self.assertEqual(
            self.example_annotation.get_barcodes5(), (["NNNNGTAACNNNNN", "NNNNGTAACNNNNN", "NNNNCCGGANNN"], 0)
        )

    def test_get_adapters(self):
        """Test if the adapter sequences are returned."""
        self.assertEqual(
            self.example_annotation.get_adapters(),
            ["AGATCGGAAGAGCGGTTCAG", "AGATCGGAAGAGCGGTTCAG", "AGATCGGAAGAGCGGTTCAG"],
        )

    def test_validate_barcode_uniqness(self):
        """Test the validation of barcode uniqueness."""
        linker_column_name = "Linker"
        barcode5_column_name = "5' barcode"

        # Test if correct content passes validation.
        self.example_annotation.validate_barcode_uniqness()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation passes when no linkers are given but the 5' barcodes are unique.
        self.change_column_content(
            column_name=barcode5_column_name,
            replacement=["NNNN,GTAAC_0,NNNNN", "NNNN,GTCAC_0,NNNNN", "NNNN,CCGGA_0,NNN"],
        )
        self.change_column_content(column_name=linker_column_name, replacement=["L3", "L3", "L3"])

        self.example_annotation.validate_barcode_uniqness()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when no linkers are given but the 5' barcodes are not unique.
        self.change_column_content(
            column_name=barcode5_column_name,
            replacement=["NNNN,GTAAC_0,NNNNN", "NNNN,GTAAC_0,NNNNN", "NNNN,CCGGA_0,NNN"],
        )
        self.change_column_content(column_name=linker_column_name, replacement=["L3", "L3", "L3"])

        self.example_annotation.validate_barcode_uniqness()

        message = "Barcodes on 5' end are not unique."
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

        # Test if the validation fails when linkers are given but the combinations are not unique.
        self.change_column_content(column_name=linker_column_name, replacement=["L3-GTC", "L3-GTC", "L3"])

        self.example_annotation.validate_barcode_uniqness()

        message = "Combination of barcodes on 3' and 5' end is not unique."
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

    def test_validate_single_adapter(self):
        """Test the validation of same adapter over all samples."""
        # Test if the validation passes when all adapters are the same.
        self.example_annotation.validate_single_adapter()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when all adapters are not same.
        self.change_column_content(
            column_name="3' adapter",
            replacement=["CAGATCGGAAG_1,AGCGGTTCAG_2", "AGATCGGAAG_1,AGCGGTTCAG_2", "AGATCGGAAG_1,AGCGGTTCAG_2"],
        )

        self.example_annotation.validate_single_adapter()

        message = "All samples should have the same adapter sequence."
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

    def test_validate_sample_names(self):
        """Test the validation of unique sample names."""
        # Test if the validation passes when all sample names are unique and present.
        self.example_annotation.validate_sample_names()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when all sample names are not unique.
        self.example_annotation.sample_names = ["Sample_1", "Sample_1", "Sample_3"]

        self.example_annotation.validate_sample_names()

        message = "Sample names should be unique, but these names are used multiple times: Sample_1."
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

        # Test if the validation fails when not all rows have sample names.
        self.example_annotation.sample_names = ["Sample_1", "", "Sample_3"]

        self.example_annotation.validate_sample_names()

        message = "One or more values for annotation field <Sample name> is missing."
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

        # Test if the validation fails when sample name is to long.
        to_long_name = "a" * (DescriptorSheet.max_sample_name_size + 1)
        self.example_annotation.sample_names = ["Sample_1", to_long_name, "Sample_3"]

        self.example_annotation.validate_sample_names()

        message = 'Sample name "{}" should be shorter than 99 characters.'.format(to_long_name)
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

    def test_validate_columns_names(self):
        """Test the validation of present columns."""
        # Test if the validation passes when all columns are present.
        self.example_annotation.validate_column_names()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when all columns are not present.
        self.example_annotation.column_names[1] = ""
        self.example_annotation.validate_column_names()

        message = "Annotation file does not contain all the required columns. Missing columns: {'Collection name'}"
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

    def test_validate_date(self):
        """Test the validation of date format."""
        date_column_name = "Date of gel images in lab notebook (optional)"

        # Test if the validation passes when the date format is correct.
        self.example_annotation.validate_date(date_column_name)

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when the date format is incorrect.
        self.change_column_content(
            column_name=date_column_name, replacement=["2020/04/10", "2020-04-10", "2020-04-10"]
        )

        self.example_annotation.validate_date(date_column_name)

        message = "SAMPLE: Sample_1 - Incorrect date format (2020/04/10), should be YYYY-MM-DD."
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

    def test_validate_yes_no(self):
        """Test the validation that ``value`` is "yes" or "no"."""
        yes_no_column_name = "Consensus mapping (optional)"

        # Test if the validation passes when value is correct.
        self.example_annotation.validate_yes_no(yes_no_column_name)

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation passes when there are no values.
        self.change_column_content(column_name=yes_no_column_name, replacement=["", "", ""])

        self.example_annotation.validate_yes_no(yes_no_column_name)

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when value is incorrect.
        self.change_column_content(column_name=yes_no_column_name, replacement=["yey", "no", ""])

        self.example_annotation.validate_yes_no(yes_no_column_name)

        message = 'SAMPLE: Sample_1 - Value yey in column {} should be "yes", "no" or empty.'.format(
            yes_no_column_name
        )
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

    def test_validate_required_columns(self):
        """Test the validation of essential input presence."""
        # Test if the validation passes when all values are given.
        self.example_annotation.validate_required_columns()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when not all values are given.
        self.change_column_content(column_name="Scientist", replacement=["Super Smarty", "", "Super Smarty"])

        self.example_annotation.validate_required_columns()

        message = "SAMPLE: Sample_2 - Value for column Scientist is missing."
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

    def test_validate_column_entry_viability(self):
        """Test that the validation recognises the entry."""
        tested_entries = "Method"

        # Test if the validation passes when the entry exists.
        self.example_annotation.validate_column_entry_viability(tested_entries, METHOD)

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when the entry does not exist.
        self.change_column_content(column_name=tested_entries, replacement=["iCLIP", "", "Magick"])

        self.example_annotation.validate_column_entry_viability(tested_entries, METHOD)

        message_one = "SAMPLE: Sample_2 -  is not a valid entry for column {}.".format(tested_entries)
        message_two = "SAMPLE: Sample_3 - Magick is not a valid entry for column {}.".format(tested_entries)
        self.assertEqual(self.example_annotation.errors, [message_one, message_two])
        self.example_annotation.errors.clear()

    def test_validate_integer_columns(self):
        """Test that validation correctly checks if ``value`` is integer or empty."""
        integer_column = "Replicate (optional)"

        # Test if the validation passes when the entries do exists and are correct.
        self.example_annotation.validate_integer_columns()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation passes when the entries do not exists.
        self.change_column_content(column_name=integer_column, replacement=["", "", ""])

        self.example_annotation.validate_integer_columns()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when the entries are not correct.
        self.change_column_content(column_name=integer_column, replacement=["one", "2", "3"])

        self.example_annotation.validate_integer_columns()

        message = "SAMPLE: Sample_1 - Value one in column {} should be integer.".format(integer_column)
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

    def test_validate_protein(self):
        """Test validation protein."""
        protein_column = "Protein"

        # Test if the validation passes when the entries are from human or mouse and are correct.
        self.example_annotation.validate_protein()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation passes when there are no entries.
        self.change_column_content(column_name=protein_column, replacement=["", "", ""])

        self.example_annotation.validate_protein()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when the entries are from human or mouse and are incorrect.
        self.change_column_content(column_name=protein_column, replacement=["no_protein", "", ""])

        self.example_annotation.validate_protein()

        message = "SAMPLE: Sample_1 - Gene symbol no_protein is either invalid or Knowledge Base cannot be reached."
        self.assertEqual(self.example_annotation.errors, [message])
        self.example_annotation.errors.clear()

    def test_validate_tissue(self):
        """Test validation of ``cells/tissue`` entries."""
        # Test if the validation passes when the entries are correct.
        self.example_annotation.validate_tissue()

        self.assertFalse(self.example_annotation.errors)
        self.example_annotation.errors.clear()

        # Test if the validation fails when the entries are not correct.
        self.change_column_content(column_name="Cells/tissue", replacement=["HEK293", "", "bone"])

        self.example_annotation.validate_tissue()

        message_one = "SAMPLE: Sample_2 -  is not a valid entry for the cells/tissue annotation field."
        message_two = "SAMPLE: Sample_3 - bone is not a valid entry for the cells/tissue annotation field."
        self.assertEqual(self.example_annotation.errors, [message_one, message_two])
        self.example_annotation.errors.clear()


if __name__ == "__main__":
    unittest.main()
