#!/bin/bash
set -e

echo " Starting installation of local tools..."

ADFR_DIR="$PWD/ADFR"

if [ -f "$ADFR_DIR/bin/prepare_receptor" ]; then
    echo "✅ ADFR Suite already installed. Skipping."
else
    echo " Downloading and installing ADFR Suite..."
    mkdir -p "$ADFR_DIR"

    wget --no-check-certificate -O ADFRsuite.tar.gz https://ccsb.scripps.edu/adfr/download/1038/
    tar -zxvf ADFRsuite.tar.gz
    cd ADFRsuite_x86_64Linux_1.0
    printf 'Y\n' | ./install.sh -d "$ADFR_DIR"
    cd ..
    rm -rf ADFRsuite.tar.gz ADFRsuite_x86_64Linux_1.0

    echo "✅ ADFR installed successfully at $ADFR_DIR"
fi

UNIDOCK_DIR="$PWD/bin"
UNIDOCK_BIN="$UNIDOCK_DIR/unidock"

if [ -f "$UNIDOCK_BIN" ]; then
    echo "✅ Uni-Dock already installed. Skipping."
else
    echo "✅ Downloading Uni-Dock executable directly..."
    mkdir -p "$UNIDOCK_DIR"

    wget --no-check-certificate -O "$UNIDOCK_BIN" "https://github.com/dptech-corp/Uni-Dock/releases/download/1.1.0/unidock-1.1.0-cuda120-linux-x86_64"
    
    chmod +x "$UNIDOCK_BIN"
    
    echo "✅ Uni-Dock installed successfully at $UNIDOCK_BIN"
fi

echo "=========================================="
echo " All tools are ready locally inside your project!"
echo " ~/.bashrc was NOT modified, so your Conda environments remain perfectly safe."