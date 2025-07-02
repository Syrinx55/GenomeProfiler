#!/usr/bin/env bash

# --- BEGIN MATERIAL UNDER LICENSE: CC BY-SA 4.0 ---
# Attribution: Dave Dopson et al.
# Source: <https://stackoverflow.com/a/246128>
# License: <https://creativecommons.org/licenses/by-sa/4.0/>
get_script_dir()
{
    local SOURCE_PATH="${BASH_SOURCE[0]}"
    local SYMLINK_DIR
    local SCRIPT_DIR
    # Resolve symlinks recursively
    while [ -L "$SOURCE_PATH" ]; do
        # Get symlink directory
        SYMLINK_DIR="$( cd -P "$( dirname "$SOURCE_PATH" )" >/dev/null 2>&1 && pwd )"
        # Resolve symlink target (relative or absolute)
        SOURCE_PATH="$(readlink "$SOURCE_PATH")"
        # Check if candidate path is relative or absolute
        if [[ $SOURCE_PATH != /* ]]; then
            # Candidate path is relative, resolve to full path
            SOURCE_PATH=$SYMLINK_DIR/$SOURCE_PATH
        fi
    done
    # Get final script directory path from fully resolved source path
    SCRIPT_DIR="$(cd -P "$( dirname "$SOURCE_PATH" )" >/dev/null 2>&1 && pwd)"
    echo "$SCRIPT_DIR"
}
# --- END MATERIAL UNDER LICENSE: CC BY-SA 4.0 ---

cd "$(get_script_dir)"
conda run -n genome-profiler python "$(get_script_dir)/genome_profiler/genome_profiler.py" --gui
