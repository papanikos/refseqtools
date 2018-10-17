import argparse
from .download_domain import download
from .catalogs import my_awesome_func


def cli():
    parser = argparse.ArgumentParser(description="Tools for handling refseq information")
    subparsers = parser.add_subparsers(title="subcommands", dest='subparser_name')

    # Add subparsers and arguments
    download_parser = subparsers.add_parser('download')

    download_parser.add_argument('-d',
                                 dest='domain',
                                 help="One of 'viral', 'bacteria', 'invertebrate', 'archaea' "
                                      "'fungi', invertebrate', 'protozoa', 'vertebrate_mammalian', "
                                      "'vertebrate_other'",
                                 required=True
                                 )
    download_parser.add_argument('-o',
                                 dest='output_dir',
                                 help='Full path to output directory',
                                 required=True
                                 )
    download_parser.add_argument('-p',
                                 dest='no_of_processes',
                                 type=int,
                                 help='Number of processes to start.',
                                 default=4,
                                 required=False)

    download_parser.set_defaults(func=download)

    catalogs_parser = subparsers.add_parser('catalogs')
    catalogs_parser.add_argument('-c',
                                 dest='catalogs_dir',
                                 required=True)
    catalogs_parser.add_argument('-t',
                                 dest='target_release',
                                 type=int,
                                 required=True)
    catalogs_parser.set_defaults(func=my_awesome_func)

    args = parser.parse_args()

    if args.subparser_name == 'download':
        args.func(args.domain, args.output_dir, args.no_of_processes)
    elif args.subparser_name == 'catalogs':
        args.func(args.catalogs_dir, args.target_release)

