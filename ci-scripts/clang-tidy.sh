#!/bin/bash
#
# This code is part of the tng binary trajectory format.
#
# Copyright (c) 2020, by the GROMACS development team.
# TNG was orginally written by Magnus Lundborg, Daniel Spångberg and
# Rossen Apostolov. The API is implemented mainly by Magnus Lundborg,
# Daniel Spångberg and Anders Gärdenäs.
#
# Please see the AUTHORS file for more information.
#
# The TNG library is free software; you can redistribute it and/or
# modify it under the terms of the Revised BSD License.
#
# To help us fund future development, we humbly ask that you cite
# the research papers on the package.
#
# Check out http://www.gromacs.org for more information.

# This script runs clang tidy checks on modified files and
# reports/applies the necessary changes.
#
# See `clang-tidy.sh -h` for a brief usage, and docs/dev-manual/code-formatting.rst
# for more details.

# Parse command-line arguments
function usage() {
    echo "usage: clang-tidy.sh [-f|--force] [--parallel=#Jobs] [--rev=REV]"
    echo "           [--tidy=(off|check)]"
    echo "           [--warnings=<file>] [<action>]"
    echo "           [-B=<builddir>]"
    echo "<action>: (check*|diff|update)[-(index|workdir*)] (*=default)"
}

action="check-workdir"
declare -a diffargs
baserev="HEAD"
force=
tidy_mode=check
warning_file=
builddir=
concurrency=2
for arg in "$@" ; do
    if [[ "$arg" == "check-index" || "$arg" == "check-workdir" || \
          "$arg" == "diff-index" || "$arg" == "diff-workdir" || \
          "$arg" == "update-index" || "$arg" == "update-workdir" ]]
    then
        action=$arg
    elif [[ "$arg" == "check" || "$arg" == "diff" || "$arg" == "update" ]] ; then
        action=$arg-workdir
    elif [[ "$action" == diff-* ]] ; then
        diffargs+=("$arg")
    elif [[ "$arg" == --tidy=* ]] ; then
        tidy_mode=${arg#--tidy=}
        if [[ "$tidy_mode" != "off" && "$tidy_mode" != "check" ]] ; then
            echo "Unknown option: $arg"
            echo
            usage
            exit 2
        fi
    elif [[ "$arg" == "-f" || "$arg" == "--force" ]] ; then
        force=1
    elif [[ "$arg" == --parallel=* ]] ; then
        concurrency=${arg#--parallel=}
    elif [[ "$arg" == --rev=* ]] ; then
        baserev=${arg#--rev=}
    elif [[ "$arg" == --warnings=* ]] ; then
        warning_file=${arg#--warnings=}
    elif [[ "$arg" == -B=* ]] ; then
        builddir=${arg#-B=}
    elif [[ "$arg" == "-h" || "$arg" == "--help" ]] ; then
        usage
        exit 0
    else
        echo "Unknown option: $arg"
        echo
        usage
        exit 2
    fi
done

# Check that format is present
if [[ "$tidy_mode" != "off" ]]
then
    if [ -z "$RUN_CLANG_TIDY" ]
    then
        RUN_CLANG_TIDY=`git config hooks.runclangtidypath`
    fi
    if [ -z "$RUN_CLANG_TIDY" ]
    then
        echo "Please set the path to run-clang-tidy using the git hook"
        echo "git config hooks.runclangtidypath /path/to/run-clang-tidy-9.py"
        echo "or by setting an environment variable, e.g."
        echo "RUN_CLANG_TIDY=/path/to/run-clang-tidy-9.py"
        exit 2
    fi
    if ! which "$RUN_CLANG_TIDY" 1>/dev/null
    then
        echo "run-clang-tidy-9.py not found: $RUN_CLANG_TIDY"
        exit 2
    fi
fi

# Switch to the root of the source tree and check the config file
srcdir=`git rev-parse --show-toplevel`
pushd $srcdir >/dev/null || exit

# Actual processing starts: create a temporary directory
tmpdir=`mktemp -d -t gmxclangtidy.XXXXXX`

# Produce a list of changed files
# Only include files that have proper filter set in .gitattributes
internal_diff_args=
if [[ $action == *-index ]]
then
    internal_diff_args="--cached"
fi
git diff-index $internal_diff_args --diff-filter=ACMR $baserev >$tmpdir/difflist
cut -f2 <$tmpdir/difflist | \
    git check-attr --stdin filter | \
    sed -e 's/.*: filter: //' | \
    paste $tmpdir/difflist - | \
    grep -E '(complete_formatting|clangformat|copyright|includesort)$' >$tmpdir/filtered
cut -f2 <$tmpdir/filtered >$tmpdir/filelist_all
grep -E '(complete_formatting|clangformat)$' <$tmpdir/filtered | \
    cut -f2 >$tmpdir/filelist_clangtidy
git diff-files --name-only | grep -Ff $tmpdir/filelist_all >$tmpdir/localmods

# Extract changed files to a temporary directory
mkdir $tmpdir/org
if [[ $action == *-index ]] ; then
    git checkout-index --prefix=$tmpdir/org/ --stdin <$tmpdir/filelist_all
else
    rsync --files-from=$tmpdir/filelist_all -a $srcdir/ $tmpdir/org/ 
fi
# check for the existence of the compile_commands.json file and abort
# if it is not present. If we don't have a build directory, try the
# current source directory.
if [ -z $builddir ] ; then
    builddir=$srcdir
fi
if [[ ! -f $builddir/compile_commands.json ]] ; then
    echo "Could not find compile_commands.json in builddir=$builddir"
    echo "Make sure you gave a correct build tree and that it contains the file!"
else
    # Need to have compilation database file available somewhere above where we are using it
    rsync -a $builddir/compile_commands.json $tmpdir/org
fi
# Prepare directory to use for comparing changed and original files
cp -r $tmpdir/org $tmpdir/new

# Create output file for what was done (in case no messages get written)
touch $tmpdir/messages

# Run clang-tidy on the temporary directory
# Can only perform clang-tidy on a non-empty list of files
cd $tmpdir/new
if [[ $tidy_mode != "off" &&  -s $tmpdir/filelist_clangtidy ]] ; then
    $RUN_CLANG_TIDY `cat $tmpdir/filelist_clangtidy` -header-filter=.* -j $concurrency -fix -quiet -extra-arg=--cuda-host-only -extra-arg=-nocudainc>$tmpdir/clang-tidy.out 2>&1
    awk '/warning/,/clang-tidy|^$/' $tmpdir/clang-tidy.out | grep -v "warning.*generated." | grep -v "Suppressed .* warnings" | grep -v "clang-analyzer"  | grep -v "to display errors from all non" | sed '/^\s*$/d' > $tmpdir/clang-tidy-warnings.out
    grep '\berror:' $tmpdir/clang-tidy.out > $tmpdir/clang-tidy-errors.out || true
    if [ -s $tmpdir/clang-tidy-errors.out ]; then
        echo "Running of clang-tidy failed. Check output below for errors:"
        cat $tmpdir/clang-tidy-errors.out
        rm -rf $tmpdir
        exit 2
    fi
    # Find the changed files if necessary
    if [[ $action != diff-* ]] ; then
        msg="found code issues"
        if [[ $action == update-* ]] ; then
            msg="clang-tidy performed"
        fi
        rsync --files-from=$tmpdir/filelist_all -a $srcdir/ ./
        rsync -a $tmpdir/org/ $srcdir/
        git diff --no-index --name-only ../org/ . | \
            awk -v msg="$msg" '{sub(/.\//,""); print $0 ": " msg}' >> $tmpdir/messages
    fi
    # TODO: Consider checking whether rerunning clang-tidy causes additional changes
fi

cd $tmpdir

# If a diff was requested, show it and we are done
if [[ $action == diff-* ]] ; then
    git diff --no-index --no-prefix "${diffargs[@]}" org/ new/
    rm -rf $tmpdir
    exit 0
fi

# Find the changed files
git diff --no-index --name-only --exit-code org/ new/ | \
    sed -e 's#new/##' > $tmpdir/changed
changes=
if [[ -s $tmpdir/changed ]]
then
    changes=1
fi

# Check if changed files have changed outside the index
if grep -Ff $tmpdir/localmods $tmpdir/changed > $tmpdir/conflicts
then
    awk '{print $0 ": has changes in work tree"}' $tmpdir/conflicts \
        >> $tmpdir/messages
    if [[ ! $force && $action == update-* ]] ; then
        echo "Modified files found in work tree, skipping update. Use -f to override."
        echo "The following would have been done:"
        sort $tmpdir/messages
        rm -rf $tmpdir
        exit 2
    fi
fi

# Update the index/work tree if requested
if [[ $action == update-index ]] ; then
    grep -Ff $tmpdir/changed $tmpdir/filtered > $tmpdir/tohash
    cd $tmpdir/new
    IFS='
'
    for change in `cut -f2 $tmpdir/tohash | \
                   git --git-dir=$srcdir/.git hash-object -w --stdin-paths --no-filters | \
                   paste - $tmpdir/tohash`
    do
        # NOTE: the patterns below contain literal tabs
        sha1=${change%%	*}
        rest=${change#*	}
        mode=${rest:8:6}
        path=${rest#*	}
        path=${path%%	*}
        # Contains a literal tab
        echo "$mode $sha1	$path" >> $tmpdir/toindex
    done
    unset IFS
    git --git-dir=$srcdir/.git update-index --index-info < $tmpdir/toindex
elif [[ $action == update-workdir ]] ; then
    rsync --files-from=$tmpdir/changed $tmpdir/new/ $srcdir/
fi

# Get back to the original directory
popd >/dev/null

# Report what was done
if [ -s $tmpdir/clang-tidy-warnings.out ] ; then
    cat $tmpdir/clang-tidy-warnings.out | tee $warning_file
fi
sort $tmpdir/messages | tee -a $warning_file

rm -rf $tmpdir
exit $changes
