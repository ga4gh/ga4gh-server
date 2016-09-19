"""
ga2sam cli
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.cli as cli
import ga4gh.cli.client as cli_client
import ga4gh.converters as converters


class Ga2SamRunner(cli_client.SearchReadsRunner):
    """
    Runner class for the ga2vcf
    """
    def __init__(self, args):
        args.readGroupIds = args.readGroupId
        super(Ga2SamRunner, self).__init__(args)
        self._outputFile = args.outputFile
        self._binaryOutput = False
        if args.outputFormat == "bam":
            self._binaryOutput = True

    def run(self):
        samConverter = converters.SamConverter(
            self._client, readGroupId=self._readGroupIds[0],
            referenceId=self._referenceId, start=self._start, end=self._end,
            outputFileName=self._outputFile, binaryOutput=self._binaryOutput)
        samConverter.convert()


def getGa2SamParser():
    parser = cli.createArgumentParser("GA4GH SAM conversion tool")
    cli_client.addClientGlobalOptions(parser)
    cli_client.addUrlArgument(parser)
    parser.add_argument(
        "readGroupId",
        help="The ReadGroup to convert to SAM/BAM format.")
    cli_client.addPageSizeArgument(parser)
    cli_client.addStartArgument(parser)
    cli_client.addEndArgument(parser)
    parser.add_argument(
        "--referenceId", default=None,
        help="The referenceId to search over")
    parser.add_argument(
        "--outputFormat", "-O", default="sam", choices=["sam", "bam"],
        help=(
            "The format for object output. Currently supported are "
            "'sam' (default), which is a text-based format and "
            "'bam', which is the binary equivalent"))
    cli.addOutputFileArgument(parser)
    return parser


def ga2sam_main():
    parser = getGa2SamParser()
    args = parser.parse_args()
    if "baseUrl" not in args:
        parser.print_help()
    else:
        runner = Ga2SamRunner(args)
        runner.run()
