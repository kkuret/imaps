"""Batch download of files from iMaps server.

Examples
========

Download all possible FASTQ files from collection "ABC" into directory
``/home/john/iclip/``::

    python batch_download.py
        --collection ABC
        --types fastq
        --directory /home/john/iclip/

Download all bedgraph files from xlsites and group objects. Also download all
gene summaries from collection "ABC"::

    python batch_download.py
        --collection ABC
        --types bedgraph-xlsites,bedgraph-group,gene_summary

Full set of options for ``--types`` argument can be viewed in
``SUPPORTED_TYPES`` list.

"""
import argparse
import os
import pathlib

import resdk

SERVER_URL = "https://imaps.genialis.com"
SUPPORTED_TYPES = [
    'all',
    'fastq',
    'bam',
    'bed',
    'bedgraph',
    'bed-annotate',
    'bedgraph-annotate',
    'bed-clusters',
    'bedgraph-clusters',
    'bed-group',
    'bedgraph-group',
    'bed-xlsites',
    'bedgraph-xlsites',
    'bed_multi',
    'peaks',
    'bedgraph-peaks',
    'type_summary',
    'subtype_summary',
    'gene_summary',
    'bed-paraclu',
    'results-rnamaps',
    'results-kmers',  # New kmers version
    'figures-kmers',  # Old kmers version
    'pos_distribution-kmers',  # Old kmers version
    'kmer_counts-kmers',  # Old kmers version
]


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--collection', required=True, help="Collection name.")
    parser.add_argument(
        '-t',
        '--types',
        required=True,
        help="Types of files to download. If multiple types are given, separate them by comma.",
    )
    parser.add_argument(
        '-d',
        '--directory',
        default=None,
        help="Directory into which to download files. If not given, download to current working "
             "directory.",
    )
    return parser.parse_args()


def parse_types(input_types):
    """Parse types argument."""
    if input_types == 'all':
        return SUPPORTED_TYPES[:]

    types = []
    for type_ in input_types.split(','):
        type_ = type_.strip()
        if type_ not in SUPPORTED_TYPES:
            raise ValueError('Type "{}" is not supported.'.format(type_))
        types.append(type_.split('-'))

    return types


def get_unexisting_name(name, directory):
    """Get unexisting name it the one given already exists."""
    extension = ''.join(pathlib.Path(name).suffixes)
    basename = os.path.basename(name)[:-len(extension)]

    i = 1
    while name in os.listdir(directory):
        name = '{} ({}){}'.format(basename, i, extension)
        i += 1
    return name


def rename_if_clashing(name, directory):
    """Rename file if it alrady exists."""
    if name in os.listdir(directory):
        os.rename(
            os.path.join(directory, name),
            os.path.join(directory, get_unexisting_name(name, directory)),
        )


def main():
    """Invoke when run directly as a program."""
    args = parse_arguments()

    res = resdk.Resolwe(url=SERVER_URL)
    res.login()
    collection = res.collection.get(name=args.collection)

    types = parse_types(args.types)
    for data in collection.data:
        if data.status != 'OK':
            continue

        for type_ in types:
            # type is a tuple of size 1 or 2: (field_name) or (field_name, process_type)
            if len(type_) == 2:
                if not data.process.type.strip(':').endswith(type_[1]):
                    continue

            field_name = type_[0]

            if field_name not in data.output:
                continue

            if isinstance(data.output[field_name], list):
                for item in data.output[field_name]:
                    # Check if file name of the file to-be-downloaded will be
                    # clashing with existing filenames in download direcory. If
                    # so, rename existing file to unexisting name.
                    original_name = os.path.basename(item['file'])
                    rename_if_clashing(original_name, args.directory)
            else:
                original_name = os.path.basename(data.output[field_name]['file'])
                rename_if_clashing(original_name, args.directory)

            print('Downloading {} output of data {} ...'.format(field_name, data.name))
            data.download(field_name=field_name, download_dir=args.directory)


if __name__ == "__main__":
    main()
