"""
Convert a BAM file to a small BAM file
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.common.utils as utils


@utils.Timed()
def main():
    tool = utils.AlignmentFileTool(
        utils.AlignmentFileConstants.BAM,
        utils.AlignmentFileConstants.BAM)
    tool.parseArgs()
    tool.convert()


if __name__ == '__main__':
    main()
