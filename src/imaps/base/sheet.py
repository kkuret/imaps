"""Parse and validate iClip annotation template."""
import collections
import csv
import gzip
import re

import resdk
import xlrd
from imaps.base.constants.assets import SPECIES
from imaps.base.constants.sheet import COLUMNS, INTEGER_COLUMNS, METHOD, PROTEIN, REQUIRED_COLUMNS, TISSUE
from imaps.base.exceptions import ValidationError
from imaps.base.validation import validate_date, validate_string


class DescriptorSheet:
    """iClip annotation template operations."""

    max_sample_name_size = 99

    def __init__(self, filename):
        """Define all instance variables."""
        self.fname = filename

        self.errors = []

        self.extension = None
        self.get_extension()

        self.annotation_tab_name = self.fname.replace(self.extension, "") + "tab.gz"

        self.content = []
        self.column_names = []
        self.parse()

        self.sample_names = self.get_column("Sample name")

    def parse(self):
        """Parse annotation content."""
        if self.extension not in ("xls", "xlsx", "tab.gz"):
            self.error("File extension not recognized.")
            self.fail()

        if self.extension in ("xls", "xlsx"):
            workbook = xlrd.open_workbook(self.fname)
            worksheet = workbook.sheets()[0]
            column_names = worksheet.row_values(0)
            content = []
            for rownum in range(1, worksheet.nrows):
                ascii_row_content = []
                for cell_content, column_name in zip(worksheet.row_values(rownum), column_names):
                    # Handle non-ascii charachters:
                    try:
                        ascii_value = bytes(str(cell_content).encode("utf-8")).decode("ascii", "strict")
                    except (UnicodeEncodeError, UnicodeDecodeError):
                        for position, char in enumerate(cell_content):
                            if ord(char) > 127:
                                break
                        self.warning(
                            'Problem decoding row {}, column "{}" at position {}.'.format(
                                rownum + 1, column_name, position
                            )
                        )
                        ascii_value = bytes(str(cell_content).encode("utf-8")).decode("ascii", "ignore")
                    finally:
                        # Strip cell values as this can cause problems downstream.
                        ascii_value.strip()
                        ascii_row_content.append(ascii_value)
                content.append(ascii_row_content)

        if self.extension == "tab.gz":
            with gzip.open(self.fname, "rt") as unzipped_annotation:
                annotation_content = list(csv.reader(unzipped_annotation, delimiter="\t"))
                column_names = annotation_content[0]
                content = annotation_content[1:]

        self.content = content
        self.column_names = column_names

    def write_tab_file(self, path=None):
        """Write a compressed tab delimited file."""
        if path is None:
            path = self.annotation_tab_name

        with gzip.open(path, "wt") as tab_file:
            tab_writer = csv.writer(tab_file, delimiter="\t")
            tab_writer.writerow(self.column_names)
            tab_writer.writerows(self.content)

    def get_extension(self):
        """Obtain the full extension."""
        extension = self.fname.split(".")[-1]

        if extension == "gz":
            extension = ".".join(self.fname.split(".")[-2:])

        self.extension = extension

    def get_column(self, column_name):
        """Obtain named column content."""
        column_index = next(i for i, name in enumerate(self.column_names) if name == column_name)
        return [row[column_index] for row in self.content]

    def get_element(self, column_name, sample_name):
        """Obtain content of specific sample in the named column."""
        if column_name not in self.column_names:
            raise ValueError(f'There is no column with column name: "{column_name}"')
        if sample_name not in self.sample_names:
            raise ValueError(f'There is no sample with sample name: "{sample_name}"')

        for row, sample in zip(self.content, self.sample_names):
            for element, column in zip(row, self.column_names):
                if sample == sample_name and column == column_name:
                    return element

    def get_barcodes3(self):
        """Get 3' barcodes and number of tolerated mismatches."""
        barcodes3 = []
        mismatches = []
        for barcode in self.get_column("Linker"):
            barcodes3.append(re.sub(r"[^ACGTN]+", "", barcode))
            mismatches.append(int(re.sub(r"[^\d]+", "", barcode)))

        assert len(set(mismatches)) == 1
        return barcodes3, mismatches[0]

    def get_barcodes5(self):
        """Get 5' barcodes and number of tolerated mismatches."""
        barcodes5 = []
        mismatches = []
        for barcode in self.get_column("5' barcode"):
            barcodes5.append(re.sub(r"[^ACGTN]+", "", barcode))
            mismatches.append(int(re.sub(r"[^\d]+", "", barcode)))
        assert len(set(mismatches)) == 1
        return barcodes5, mismatches[0]

    def get_adapters(self):
        """Get adapter sequences."""
        adapters = []
        for adapter in self.get_column("3' adapter"):
            adapters.append(re.sub(r"[^ACGTN]+", "", adapter))
        return adapters

    def validate(self):
        """Run all validation functions."""
        self.validate_column_names()
        self.validate_required_columns()
        self.validate_integer_columns()

        self.validate_sample_names()
        self.validate_column_entry_viability("Method", METHOD)
        self.validate_protein()
        self.validate_tissue()
        self.validate_column_entry_viability("Species", SPECIES)
        self.validate_barcode_uniqness()
        self.validate_single_adapter()
        self.validate_yes_no("Consensus mapping (optional)")
        self.validate_date("Date of gel images in lab notebook (optional)")

        if self.errors:
            self.fail()

    def validate_barcode_uniqness(self):
        """Validate uniqness of barcodes."""
        barcodes3, _ = self.get_barcodes3()
        barcodes5, _ = self.get_barcodes5()
        if all([brc == "" for brc in barcodes3]):
            # No 3' barcodes are given, check only for uniqness of 5' ones.
            if len(barcodes5) != len(set(barcodes5)):
                self.error("Barcodes on 5' end are not unique.")
        else:
            combined = list(zip(barcodes5, barcodes3))
            if len(combined) != len(set(combined)):
                self.error("Combination of barcodes on 3' and 5' end is not unique.")

    def validate_single_adapter(self):
        """Validate that all samples have same adapter."""
        adapters = self.get_adapters()
        if len(set(adapters)) > 1:
            self.error("All samples should have the same adapter sequence.")

    def validate_sample_names(self):
        """Validate that all samples have names and are unique."""
        if len(self.sample_names) > len(set(self.sample_names)):
            repeated = [name for name, count in collections.Counter(self.sample_names).items() if count >= 2]
            repeated = ", ".join(repeated)
            self.error("Sample names should be unique, but these names are used multiple times: {}.".format(repeated))

        for name in self.sample_names:
            if not name:
                self.error("One or more values for annotation field <Sample name> is missing.")
            if len(name) >= self.max_sample_name_size:
                self.error(
                    'Sample name "{}" should be shorter than {} characters.'.format(name, self.max_sample_name_size)
                )

    def validate_column_names(self):
        """Validate if all columns are present."""
        missing_columns = COLUMNS - set(self.column_names)
        if missing_columns:
            self.error(
                "Annotation file does not contain all the required columns. Missing columns: {}".format(
                    missing_columns
                )
            )

    def validate_date(self, column_name):
        """Validate date format."""
        for sample_name, date in zip(self.sample_names, self.get_column(column_name)):
            try:
                validate_date(date, allow_empty=True)
            except ValueError:
                self.error("SAMPLE: {} - Incorrect date format ({}), should be YYYY-MM-DD.".format(sample_name, date))

    def validate_yes_no(self, column_name):
        """Validate that ``value`` is "yes" or "no"."""
        for sample_name, element in zip(self.sample_names, self.get_column(column_name)):
            try:
                validate_string(element, choices=["yes", "no"], allow_empty=True)
            except ValueError:
                self.error(
                    'SAMPLE: {} - Value {} in column {} should be "yes", "no" or empty.'.format(
                        sample_name, element, column_name
                    )
                )

    def validate_required_columns(self):
        """Check if essential input is given."""
        for column_name in REQUIRED_COLUMNS:
            for sample_name in self.sample_names:
                element = self.get_element(column_name=column_name, sample_name=sample_name)

                if not element:
                    self.error("SAMPLE: {} - Value for column {} is missing.".format(sample_name, column_name))

    def validate_integer_columns(self):
        """Validate that ``value`` is integer or empty."""
        for column_name in INTEGER_COLUMNS:
            for sample_name in self.sample_names:
                element = self.get_element(column_name=column_name, sample_name=sample_name)

                if element:
                    try:
                        int(float(element))
                    except ValueError:
                        self.error(
                            "SAMPLE: {} - Value {} in column {} should be integer.".format(
                                sample_name, element, column_name
                            )
                        )

    def validate_tissue(self):
        """Validate that ``cells/tissue`` is existent."""
        for sample_name in self.sample_names:
            tissue = self.get_element(column_name="Cells/tissue", sample_name=sample_name)

            if self.get_part_before_colon(tissue) not in TISSUE:
                self.error(
                    "SAMPLE: {} - {} is not a valid entry for the cells/tissue annotation field.".format(
                        sample_name, tissue
                    )
                )

    def validate_column_entry_viability(self, column, possible):
        """Validate that the entry is recognized."""
        for sample_name in self.sample_names:
            element = self.get_element(column_name=column, sample_name=sample_name)

            if element not in possible:
                self.error("SAMPLE: {} - {} is not a valid entry for column {}.".format(sample_name, element, column))

    def validate_protein(self):
        """Only validate protein names if species is human or mouse."""
        res = resdk.Resolwe(url="https://app.genialis.com")

        for sample_name in self.sample_names:
            species = self.get_element(column_name="Species", sample_name=sample_name)
            protein = self.get_element(column_name="Protein", sample_name=sample_name)
            gene_symbol = self.get_part_before_colon_hypen(protein)

            if gene_symbol and gene_symbol not in PROTEIN:
                if species in ["Homo sapiens", "Mus musculus"]:
                    kb_gene = res.feature.filter(source="UCSC", feature_id=[gene_symbol])
                    if not kb_gene:
                        self.error(
                            "SAMPLE: {} - Gene symbol {} is either invalid or "
                            "Knowledge Base cannot be reached.".format(sample_name, protein)
                        )

    def error(self, string):
        """Save error messages."""
        self.errors.append(string)

    def fail(self):
        """Report collected error messages."""
        raise ValidationError(self.errors)

    @staticmethod
    def get_part_before_colon(string):
        """Return part of string before first first colon."""
        try:
            return re.match(r".*?([^:]+)[:]?", string).group(1)
        except AttributeError:
            return string

    @staticmethod
    def get_part_before_colon_hypen(string):
        """Return part of string before first first colon / hypen."""
        try:
            return re.match(r".*?([^:\-]+)[:\-]?", string).group(1)
        except AttributeError:
            return string

    @staticmethod
    def warning(string):
        """Return an warning message."""
        print(string)
