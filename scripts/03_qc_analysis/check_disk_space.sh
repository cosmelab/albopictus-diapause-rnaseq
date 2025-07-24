#!/bin/bash

echo "=== Disk Space Overview ==="
df -h

echo -e "\n=== Home Directory Usage ==="
du -sh ~/.* 2>/dev/null | sort -hr | head -n 10

echo -e "\n=== Large Directories in Home ==="
du -sh ~/* 2>/dev/null | sort -hr | head -n 10

echo -e "\n=== Large Files in Home (over 1GB) ==="
find ~ -type f -size +1G -exec du -h {} \; 2>/dev/null | sort -hr

echo -e "\n=== Docker Space Usage ==="
docker system df -v 2>/dev/null || echo "Docker not installed or not running"

echo -e "\n=== Downloads Directory ==="
du -sh ~/Downloads 2>/dev/null

echo -e "\n=== Cache Directories ==="
du -sh ~/Library/Caches 2>/dev/null
du -sh ~/Library/Application\ Support 2>/dev/null

echo -e "\n=== Temporary Files ==="
du -sh /tmp 2>/dev/null

echo -e "\n=== Recommendations ==="
echo "1. Check ~/Downloads for large files"
echo "2. Clear browser caches"
echo "3. Remove old Docker images and containers"
echo "4. Clear system caches"
echo "5. Remove old application data"
echo "6. Check for large log files"
echo "7. Remove old Time Machine backups"
echo "8. Clear mail attachments and downloads" 