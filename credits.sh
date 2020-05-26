#!/bin/bash

cd tests

python -m duecredit test_import.py

printf "\n\nBIBTEX\n------\n\n"

duecredit summary --format bibtex