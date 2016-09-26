"""
ga2vcf cli
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.cli as cli
import ga4gh.cli.client as cli_client
import ga4gh.converters as converters


class Ga2VcfRunner(cli_client.SearchVariantsRunner):
    """
    Runner class for the ga2vcf
    """
    def __init__(self, args):
        super(Ga2VcfRunner, self).__init__(args)
        self._outputFile = args.outputFile
        self._binaryOutput = False
        if args.outputFormat == "bcf":
            self._binaryOutput = True

    def run(self):
        variantSet = self._client.get_variant_set(self._variantSetId)
        iterator = self._client.search_variants(
            start=self._start, end=self._end,
            reference_name=self._referenceName,
            variant_set_id=self._variantSetId,
            call_set_ids=self._callSetIds)
        # do conversion
        vcfConverter = converters.VcfConverter(
            variantSet, iterator, self._outputFile, self._binaryOutput)
        vcfConverter.convert()


def getGa2VcfParser():
    parser = cli.createArgumentParser((
        "GA4GH VCF conversion tool. Converts variant information "
        "stored in a GA4GH repository into VCF format."))
    cli_client.addClientGlobalOptions(parser)
    cli.addOutputFileArgument(parser)
    cli_client.addUrlArgument(parser)
    parser.add_argument("variantSetId", help="The variant set to convert")
    parser.add_argument(
        "--outputFormat", "-O", choices=['vcf', 'bcf'], default="vcf",
        help=(
            "The format for object output. Currently supported are "
            "'vcf' (default), which is a text-based format and "
            "'bcf', which is the binary equivalent"))
    cli_client.addReferenceNameArgument(parser)
    cli_client.addCallSetIdsArgument(parser)
    cli_client.addStartArgument(parser)
    cli_client.addEndArgument(parser)
    cli_client.addPageSizeArgument(parser)
    return parser


def ga2vcf_main():
    parser = getGa2VcfParser()
    args = parser.parse_args()
    if "baseUrl" not in args:
        parser.print_help()
    else:
        runner = Ga2VcfRunner(args)
        runner.run()
