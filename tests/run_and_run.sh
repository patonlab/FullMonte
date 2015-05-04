#!/usr/bin/env bash
while :
do
	clear
	echo "Running tests..."
	python test.py
	echo -e "\n\nPress q to quit or any other key to run again..." 
	read -s -n1 choice
	case $choice in
		q)
			echo -e "\nQuitting"
			exit 0
			;;
	esac
done