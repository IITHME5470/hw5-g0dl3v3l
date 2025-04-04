#!/bin/bash

REMOTE_USER="ee21btech11015"
REMOTE_IP="192.168.172.48"
REMOTE_DIR="/clhome/ee21btech11015/hw5"

LOCAL_PARENT="/Users/dheerajmahendiran/Desktop/parallelComputing/hw5-g0dl3v3l"

echo "📥 Pulling directory from $REMOTE_USER@$REMOTE_IP:$REMOTE_DIR to $LOCAL_PARENT"

# This will create hw5 inside hw5-g0dl3v3l
scp -r \
  HostKeyAlgorithms=ssh-rsa \
  PubkeyAcceptedAlgorithms=ssh-rsa \
  "$REMOTE_USER@$REMOTE_IP:$REMOTE_DIR" "$LOCAL_PARENT"

if [ $? -eq 0 ]; then
  echo "✅ Transfer successful!"
else
  echo "❌ Transfer failed!"
fi
