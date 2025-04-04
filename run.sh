#!/bin/bash

# === CONFIGURATION ===

# Local directory you want to send
LOCAL_DIR="/Users/dheerajmahendiran/Desktop/parallelComputing/hw5-g0dl3v3l/hw5"

# Remote user and IP
REMOTE_USER="ee21btech11015"
REMOTE_IP="192.168.172.48"

# Remote target directory (destination on server)
REMOTE_DIR="/clhome/ee21btech11015"

# === TRANSFER ===

echo "Transferring directory $LOCAL_DIR to $REMOTE_USER@$REMOTE_IP:$REMOTE_DIR"

scp -r -o HostKeyAlgorithms=+ssh-rsa -o PubkeyAcceptedAlgorithms=+ssh-rsa "$LOCAL_DIR" \
    "${REMOTE_USER}@${REMOTE_IP}:${REMOTE_DIR}"

# === DONE ===

if [ $? -eq 0 ]; then
  echo "✅ Transfer successful!"
else
  echo "❌ Transfer failed!"
fi
