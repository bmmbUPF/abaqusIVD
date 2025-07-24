# Open the given file for reading and a new file for writing
with open('TransportDisc_elset_NP_Z5.inp', 'r') as infile, open('output_file.txt', 'w') as outfile:
    # Read each line from the input file
    for line in infile:
        # Split the line based on commas to get individual numbers
        numbers = line.strip().split(',')
        # Increment each number by 30,000
        incremented_numbers = [str(int(num.strip()) + 30000) for num in numbers]
        # Join the incremented numbers back with commas and write to the output file
        outfile.write(', '.join(incremented_numbers) + '\n')

print("Updated file has been written to output_file.txt")

