#!/bin/bash

#   Check changed vignettes for SeuratWrappers
#   Set global options to cause the script to fail upon failure of any one step
set -eo pipefail

#   A simple function to get the extension of files
function extension() {
    local fname="${1}" # Name of file to get extension of
    echo "$(echo ${fname} | rev | cut -f 1 -d '.' | rev)"
}

export -f extension

#   Get the remote for Satija Lab
SATIJA_BRANCH="$(Rscript -e "cat(Seurat:::RandomName())")"
(set -x; git fetch https://github.com/satijalab/seurat-wrappers +master:"${SATIJA_BRANCH}")

#   Get the branch for this PR
PR_HASH="$(set -x; git log -n1 --format=format:'%H')"

#   Get differences between files
declare -a DIFFS=($(set -x; git diff --name-only "${PR_HASH}" "${SATIJA_BRANCH}"))

#   Figure out which files have corresponding Rmds
declare -a DIFF_RMDS=()
declare -a MISSING=()
for DFILE in ${DIFFS[@]}; do
    case $(dirname ${DFILE}) in # Only certain files will be checked
        docs) # Ensure we're only checking changed vignettes; these will have an extension of rmd or Rmd
            $(echo $(extension ${DFILE}) | grep -iw rmd > /dev/null 2> /dev/null) || continue
            DIFF_RMDS+=("${DFILE}")
            ;;
        R) # If a source file has changed, ensure it corresponds to a vignette
            if [[ $(basename ${DFILE}) == 'internal.R' ]]; then
                #   If internal.R has changed, throw a warning for manual checks, but don't check everything on Azure
                echo "WARNING: internal.R has changed, please check all vignettes" >&2
                continue
            elif [[ $(echo $(extension ${DFILE}) | grep -iw r > /dev/null 2> /dev/null; echo "$?") -eq 0 ]]; then
                #   If an R file has changed, ask if it has a vignette
                #   If so, add the vignette to the list of vignettes to check
                #   If not, put it in list of source files missing a vignette
                BNAME="$(basename ${DFILE} .$(extension ${DFILE}))"
                DVIGNETTE=$(find docs -maxdepth 1 -iregex "^docs/${BNAME}\.rmd")
                [[ "${#DVIGNETTE}" -eq 0 ]] && MISSING+=("${DFILE}") || DIFF_RMDS+=("${DVIGNETTE}")
            else
                #   Non-R files shouldn't be here
                continue
            fi
            ;;
        *) # All other files are not checked
            continue
            ;;
    esac
done

#   Do we have vignettes to check
if [[ ${#DIFF_RMDS[@]} -eq 0 && ${#MISSING[@]} -gt 0 ]]; then
    #   No, but source files changed, throw an error
    echo "ERROR: Missing vignettes for all changed source files" >&2
    exit 1
elif [[ ${#MISSING[@]} -gt 0 ]]; then
    #   Yes, but some source files changed and we couldn't find a vignette
    echo -e "WARNING: Missing vignettes for the following source files:" >&2
    for MV in ${MISSING[@]}; do echo -e "\t${MV}" >&2; done
elif [[ ${#DIFF_RMDS[@]} -eq 0 ]]; then
    #   No, and no source files changed
    echo "No changed vignettes" >&2
    exit 0
else
    :
fi

#   Store new vignettes in test-build dir
mkdir test-build

#   Filter our changed vignettes list to only unique vignettes
declare -a UNIQ_DIFFS=($(echo ${DIFF_RMDS[@]} | tr ' ' '\n' | sort | uniq))
for I in $(seq 1 ${#UNIQ_DIFFS[@]}); do
    TFILE="${UNIQ_DIFFS[$((${I} - 1))]}"
    echo "Testing vignette ${TFILE} (vignette ${I} of ${#UNIQ_DIFFS[@]})" >&2
    (set -x; Rscript -e "rmarkdown::render('${TFILE}', output_format = 'all', output_dir = 'test-build')")
done
