#!/usr/bin/env python

import numpy as np
import itertools
import sys
import tracemalloc

sys.path.append('../')
from fock_basis import get_fock_full_N, get_fock_full_N_dev

def print_memory_usage(func):
    """
    Decorator to print memory usage for each function.
    """
    def wrapper(*args, **kwargs):
        # Take a snapshot before running the function
        snapshot_before = tracemalloc.take_snapshot()

        # Run the function
        result = func(*args, **kwargs)

        # Take a snapshot after running the function
        snapshot_after = tracemalloc.take_snapshot()

        # Calculate the difference and print memory usage
        stats = snapshot_after.compare_to(snapshot_before, "lineno")
        print(f"\n[ Memory Usage for {func.__name__} ]")
        for stat in stats[:10]:  # print top 10 memory changes
            print(stat)

        return result

    return wrapper

# Start tracing memory allocation
tracemalloc.start()

# Set parameters
norb = 30
N = 4

# Measure memory usage for get_fock_full_N
snapshot_before_a = tracemalloc.take_snapshot()
a = get_fock_full_N(norb, N)
snapshot_after_a = tracemalloc.take_snapshot()

# Calculate total memory used by get_fock_full_N
stats_a = snapshot_after_a.compare_to(snapshot_before_a, "lineno")
total_memory_a = sum(stat.size_diff for stat in stats_a)

print("\n[ Memory Usage for Original function ]")
for stat in stats_a[:10]:  # print top 10 memory changes
    print(stat)
print(f"Total memory usage for get_fock_full_N: {total_memory_a} bytes")

# Measure memory usage for get_fock_full_N_dev
snapshot_before_b = tracemalloc.take_snapshot()
b = get_fock_full_N_dev(norb, N)
snapshot_after_b = tracemalloc.take_snapshot()

# Calculate total memory used by get_fock_full_N_dev
stats_b = snapshot_after_b.compare_to(snapshot_before_b, "lineno")
total_memory_b = sum(stat.size_diff for stat in stats_b)

print("\n[ Memory Usage for Dev. function ]")
for stat in stats_b[:10]:  # print top 10 memory changes
    print(stat)
print(f"Total memory usage for get_fock_full_N_dev: {total_memory_b} bytes")

# Compare total memory usage
print("\n[ Memory Usage Comparison ]")
print(f"Difference in memory usage: {total_memory_b - total_memory_a} bytes")

# Stop tracing memory
tracemalloc.stop()



