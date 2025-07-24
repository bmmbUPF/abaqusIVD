import sys

def read_elements(file_path):
    """
    Read elements from the file and return a dictionary.
    """
    elements = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            element_id = int(lines[i].split(',')[0].strip())
            node_ids = [int(node.strip()) for node in lines[i].split(',')[1:] + lines[i+1].split(',') if node.strip()]
            elements[element_id] = node_ids
    return elements

def write_elements(file_path, elements):
    """
    Write elements to a file in the specified format.
    """
    with open(file_path, 'w') as f:
        for element_id, node_ids in elements.items():
            # First line: element_id and the first 15 nodes
            f.write("{:>5}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6}, {:>6},\n".format(
                element_id, *node_ids[:16]))

            # Second line: Nodes from the 16th up to the last but one
            remaining_nodes = node_ids[15:-1]
            format_str = "       " + ", ".join(["{:>6}"] * len(remaining_nodes))  # 7 spaces before
            if remaining_nodes:
                f.write(format_str.format(*remaining_nodes))
                
            # Still on the second line: Write the last node without a trailing comma
            f.write(", {:>6}\n".format(node_ids[-1]))
                

def main():
    # Check for the correct number of arguments
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <path_to_inp_file>")
        return

    # Read elements from the input file
    file_path = sys.argv[1]
    elements = read_elements(file_path)
    
    # Find the maximum element id
    max_element_id = max(elements.keys())
    ElemOffset=30000
    print(f"Number of elements: {len(elements)}")
    print(f"Maximum element id: {max_element_id}")
    print(f"ElemOffset: {ElemOffset}")
    
    # Create new element ids
    new_elements = {old_id + ElemOffset: nodes for old_id, nodes in elements.items()}

    # Write the correspondence file
    with open('correspondence.txt', 'w') as f:
        for old_id, new_id in zip(elements.keys(), new_elements.keys()):
            f.write(f"{old_id}, {new_id}\n")

    # Write the new elements to a file
    write_elements('New_elems.inp', new_elements)
    print(f"New element ids written to TransportDisc_elements_C3D20P_new.inp")
    print(f"Element id correspondence written to correspondence.txt")

if __name__ == "__main__":
    main()
