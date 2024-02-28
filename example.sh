#!/bin/bash

./create_data -n 10k -d 4 -o points.10k.4d.bin
./create_data -n 20k -d 4 -o points.20k.4d.bin
./create_data -n 40k -d 4 -o points.40k.4d.bin
./create_data -n 80k -d 4 -o points.80k.4d.bin

./test_construction points.10k.4d.bin
./test_construction points.20k.4d.bin
./test_construction points.40k.4d.bin
./test_construction points.80k.4d.bin

./test_all points.10k.4d.bin 0.5
./test_all points.20k.4d.bin 0.5
./test_all points.40k.4d.bin 0.5
./test_all points.80k.4d.bin 0.5

