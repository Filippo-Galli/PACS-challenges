#!/bin/bash

function report_generation() {
	local output=$2

	# Extract iteration count
	iter_count=$(echo "$output" | grep -oP 'Iter: \K\d+')

	# Extract time
	total_time=$(echo "$output" | grep -oP 'Time: \K[^ ]+(?: ms)?')
	# Remove ' ms'
	total_time=${total_time% ms}

	# Extract error to exact solution
	error_to_exact_solution=$(echo "$output" | grep -oP '(?<=Error with exact solution: )[\d.]+([eE][-+]?\d+)?')

	printf "\t size: %3s Iteration count: %6d \t total_time %8d ms \t error: %8s\n" "$1" "$iter_count" "$total_time" "$error_to_exact_solution" >>report.txt
}

# Remove the old report file
rm -f report.txt

# Compile the code in test mode
make clean
make test_compile

printf "\n##############################################################################################################\n"
printf "This test are done with 8*pi^2*sin(2*pi*x)*sin(2*pi*y) as f and the correct solution u is sin(2*pi*x)*sin(2*pi*y)\n"
printf "Firstly, it will ask for the max power of 2 used as dimension of matrix to test (from 2 to 9) and then the max number of process of MPI to test.\n"
printf "All the result will be available on report.txt inside this folder.\n"
printf "##############################################################################################################\n"

# Prompt the user to enter a number for max power of 2
printf "\n1) Insert max power of 2 used as dimension of matrix to test (from 2 to 9):"
read num
# Prompt the user to enter a number for max number of MPI processes
printf "2) Insert max number of process of MPI to test (at least 2):"
read thread

# Check if the input is a valid number and greater than 1
if ! [[ "$thread" =~ ^[0-9]+$ ]] || ((thread <= 1)); then
	echo "Error: Not a thread number < 1 or not valid" >&2
	exit 1
fi

# Check if the input is a valid number and greater than 2
if ! [[ "$num" =~ ^[0-9]+$ ]] || ((num <= 2)); then
	echo "Error: Not a exponent < 2 or not valid" >&2
	exit 1
fi

# Sequential test
printf "\nStarting sequential test:\n\n"
printf "\n Processes: 1 \n" >>report.txt
for ((i = 2; i <= num; i++)); do
	size=$((2 ** i))
	printf "########### Test %d: with %d grid points ########### \n" "$((i - 2))" "$size"
	output=$(mpiexec -n 1 ./main $size "8*pi^2*sin(2*pi*x)*sin(2*pi*y)" 1 "0")

	# Extract information from the output
	report_generation $size "$output"

	printf "\n"
done

printf "\n" >>report.txt

# Parallel test
printf "\nStarting parallel test:\n\n"
for ((j = 2; j <= thread; j++)); do
	printf "\n Processes: $j \n" >>report.txt
	for ((i = 2; i <= num; i++)); do
		size=$((2 ** i))
		printf "########### Processes %d: with %d grid points ###########\n" "$j" "$size"
		output=$(mpiexec -n $j ./main $size "8*pi^2*sin(2*pi*x)*sin(2*pi*y)" 1 "0")

		# Extract information from the output
		report_generation $size "$output"

		printf "\n"
	done
	printf "\n" >>report.txt
done

