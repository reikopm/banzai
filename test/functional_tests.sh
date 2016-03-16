#!/bin/bash

TEST_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR=$(dirname $TEST_DIR)
cd $TEST_DIR

# Reiko's 18S
bash ../BOG18S_banzai.sh ./BOG_kelp18S/BOG_kelp18S_params.sh

# Kevan's 18S
