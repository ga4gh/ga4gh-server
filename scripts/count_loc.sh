#!/bin/sh
#
# Prints the lines of code in the server repository
#

# lines changed according to git
# use the --stat flag instead of --shortstat for a breakdown by file
INITIAL_COMMIT_SHA=8f32b89ddf7cbb39a0c03b8f50cae9a0d77d472c
git diff --shortstat $INITIAL_COMMIT_SHA
echo

# more comprehensive breakdown by line type and file language
# (doesn't count some things, e.g. doc files)
GIT_FILE_LIST=$(git ls-files)
cloc $GIT_FILE_LIST
