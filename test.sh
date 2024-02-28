#!/bin/bash

./create_data -n 1m -d 8 -o points.1m.8d.bin
./create_data -n 2m -d 8 -o points.2m.8d.bin
./create_data -n 4m -d 8 -o points.4m.8d.bin
./create_data -n 8m -d 8 -o points.8m.8d.bin

./test_construction points.1m.8d.bin
./test_construction points.2m.8d.bin
./test_construction points.4m.8d.bin
./test_construction points.8m.8d.bin

