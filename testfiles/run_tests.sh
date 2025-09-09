#!/bin/bash

# This script runs a series of tests for ../olanet.py

# Toy data
echo "--- Running Test 1: Identical binary trees ---"
python ../olanet.py --trees test_trees_binary_identical.nwk --dates generated_dates.csv
echo ""

echo "--- Running Test 2: Different binary trees ---"
python ../olanet.py --trees test_trees_binary_different.nwk --dates generated_dates.csv
echo ""

echo "--- Running Test 3: Resolvable multifurcations ---"
python ../olanet.py --trees test_trees_multifurcation.nwk --dates generated_dates.csv
echo ""

echo "--- Running Test 4: Random permutations ---"
python ../olanet.py --trees test_trees_multifurcation.nwk --permutations 100
echo ""

echo "--- Running Test 5: Random permutations on non-identical ---"
python ../olanet.py --trees test_trees_binary_different.nwk --permutations 100
echo ""

# Lamprologini dataset
echo "--- Running Lamprologini example (4 reticulations expected) ---"
python ../olanet.py -t Lamprologini.tre -p 1000
echo ""

echo "--- All tests completed. ---"
