#!/bin/bash

diff -s <(./test_cover_tree -s 10 -n 1000 -d 2 -r 0.25) <(./test_brute_force -s 10 -n 1000 -d 2 -r 0.25)
diff -s <(./test_cover_tree -s 20 -n 1000 -d 4 -r 0.35) <(./test_brute_force -s 20 -n 1000 -d 4 -r 0.35)
diff -s <(./test_cover_tree -s 30 -n 2000 -d 16 -r 0.6) <(./test_brute_force -s 30 -n 2000 -d 16 -r 0.6)
