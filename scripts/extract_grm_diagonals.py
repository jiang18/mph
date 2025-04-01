#!/usr/bin/env python3

import struct
import argparse
import os

def extract_diagonal_elements(input_file, output_file):
    """
    Extract diagonal elements from a binary file containing:
    1. First element: 4-byte integer equal to the number of rows of the square matrix
    2. Second element: 4-byte float weight value
    3. Remaining elements: lower-triangle matrix elements in column-major format
    """
    with open(input_file, 'rb') as f:
        # Read the first element (4-byte integer): number of rows
        matrix_size_bytes = f.read(4)
        if not matrix_size_bytes or len(matrix_size_bytes) < 4:
            raise ValueError("File too small to contain matrix size")
        matrix_size = struct.unpack('i', matrix_size_bytes)[0]
        
        # Read the second element (4-byte float): weight value
        weight_bytes = f.read(4)
        if not weight_bytes or len(weight_bytes) < 4:
            raise ValueError("File too small to contain weight value")
        weight = struct.unpack('f', weight_bytes)[0]
        
        print(f"Matrix size: {matrix_size}x{matrix_size}")
        print(f"Weight value: {weight}")
        
        # Pre-calculate all diagonal positions
        print("Calculating diagonal positions...")
        diagonal_positions = {}
        position = 8  # Start position after matrix_size and weight (8 bytes total)
        
        for col in range(matrix_size):
            # In column-major format for lower triangular matrix:
            # - We only store elements where row >= col
            # - The diagonal element (where row == col) is the first element stored for each column
            diagonal_positions[col] = position
            
            # Move to the next column
            # Skip all elements in this column (including the diagonal we just recorded)
            elements_in_column = matrix_size - col
            position += elements_in_column * 4
        
        # Verify expected file size
        expected_elements = 1 + 1 + (matrix_size * (matrix_size + 1)) // 2  # size + weight + matrix elements
        expected_file_size = expected_elements * 4  # each element is 4 bytes
        actual_file_size = os.path.getsize(input_file)
        
        if expected_file_size != actual_file_size:
            print(f"Warning: File size mismatch.")
            print(f"Expected {expected_file_size} bytes for matrix size {matrix_size}")
            print(f"Actual file size: {actual_file_size} bytes")
            print(f"Difference: {abs(expected_file_size - actual_file_size)} bytes")
        
        # Extract diagonal elements
        print("Starting extraction of diagonal elements...")
        with open(output_file, 'w') as f_out:
            # Write the diagonal values divided by weight
            
            # Read and write each diagonal element
            for col in range(matrix_size):
                # Look up the position from our hash map
                position = diagonal_positions[col]
                
                # Jump to the position of this diagonal element
                f.seek(position)
                
                # Read the 4-byte float
                element_bytes = f.read(4)
                if not element_bytes or len(element_bytes) < 4:
                    print(f"Warning: Unexpected end of file at position {position}")
                    break
                    
                element = struct.unpack('f', element_bytes)[0]
                
                # Divide by weight and write only the value
                normalized_element = element / weight
                f_out.write(f"{normalized_element}\n")
                
                # Progress indication for large matrices
                if col % 10000 == 0 and col > 0:
                    print(f"Processed {col}/{matrix_size} diagonal elements...")

def main():
    parser = argparse.ArgumentParser(description='Extract diagonal elements from MPH binary GRM file containing lower triangular in column-major format.')
    parser.add_argument('input_file', help='Path to the binary input file')
    parser.add_argument('output_file', help='Path to the text output file')
    args = parser.parse_args()
    
    try:
        print(f"Processing file: {args.input_file}")
        extract_diagonal_elements(args.input_file, args.output_file)
        print(f"Successfully extracted diagonal elements to {args.output_file}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()