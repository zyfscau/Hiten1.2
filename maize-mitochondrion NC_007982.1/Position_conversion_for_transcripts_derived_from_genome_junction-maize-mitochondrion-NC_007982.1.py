import csv
import os

# Define input and output filenames
input_file = "Filter-the-results+FL.csv"
output_file1 = "High_confidence_transcripts.csv"
output_file2 = "junction-refinement.csv"
output_file3 = "junction-refinement-result.csv"

# Check if input file exists
if not os.path.exists(input_file):
    print(f"Error: Input file '{input_file}' does not exist!")
    exit(1)

# Step 1: Read input file and classify data
selected_rows = []   # Rows that meet conditions
unselected_rows = [] # Rows that don't meet conditions

with open(input_file, 'r') as f_in:
    reader = csv.reader(f_in)
    header = next(reader)  # Read header row
    
    for row in reader:
        if len(row) < 77:  # Ensure row has enough columns
            continue
            
        col_5end = row[73]  # Column 74 (index 73)
        col_3end = row[74]  # Column 75 (index 74)
        try:
            op_val = float(row[75])  # Column 76 (index 75)
        except ValueError:
            op_val = 0.0
        col_strand = row[76]  # Column 77 (index 76)
        
        # Check four conditions
        condition1 = (col_5end == "569630") and (op_val < 0) and (col_strand == "-")
        condition2 = (col_5end == "1") and (op_val < 0) and (col_strand == "+")
        condition3 = (col_3end == "1") and (op_val < 0) and (col_strand == "-")
        condition4 = (col_3end == "569630") and (op_val < 0) and (col_strand == "+")
        
        if condition1 or condition2 or condition3 or condition4:
            selected_rows.append(row.copy())  # Use copy
        else:
            unselected_rows.append(row)

# Step 2: Write output files
# Write unselected rows to output_file1
with open(output_file1, 'w', newline='') as f_out1:
    writer = csv.writer(f_out1)
    writer.writerow(header)
    writer.writerows(unselected_rows)

# Write selected rows to output_file2 (original data)
with open(output_file2, 'w', newline='') as f_out2:
    writer = csv.writer(f_out2)
    writer.writerow(header)
    writer.writerows(selected_rows)

# Step 3: Process selected rows - fix numerical output to integer format
modified_rows = []
for row in selected_rows:
    # Copy row for modification
    new_row = row.copy()
    
    col_5end = new_row[73]
    col_3end = new_row[74]
    try:
        op_val = float(new_row[75])
    except ValueError:
        op_val = 0.0
    col_strand = new_row[76]
    
    # Apply modification rules - ensure output is integer
    if col_5end == "569630" and op_val < 0 and col_strand == "-":
        # Output absolute value of OP to 5'end (integer format)
        new_row[73] = str(int(abs(op_val)))  # Convert to integer
    elif col_3end == "569630" and op_val < 0 and col_strand == "+":
        # Output absolute value of OP to 3'end (integer format)
        new_row[74] = str(int(abs(op_val)))  # Convert to integer
    elif col_5end == "1" and op_val < 0 and col_strand == "+":
        # Calculate new value and convert to integer
        new_value = 1 + 569630 + op_val
        new_row[73] = str(int(new_value))
    elif col_3end == "1" and op_val < 0 and col_strand == "-":
        # Calculate new value and convert to integer
        new_value = 1 + 569630 + op_val
        new_row[74] = str(int(new_value))
    
    # Set OP to 0 (integer format)
    new_row[75] = "0"
    
    # Create connection string (column 78) - using integer format values
    connected_str = f"MT|{new_row[73]}|{new_row[74]}|{new_row[75]}|{col_strand}"
    
    # Add to row as new column
    new_row.append(connected_str)
    modified_rows.append(new_row)

# Step 4: Write modified data to output_file3 (including original columns + new column)
new_header = header + ["Connected"]
with open(output_file3, 'w', newline='') as f_out3:
    writer = csv.writer(f_out3)
    writer.writerow(new_header)
    writer.writerows(modified_rows)

# Step 5: Perform replacement operation in first column
# Create list containing only replaced data (for final output)
final_output_rows = []

# Process header row (keep as is)
final_output_rows.append(header.copy())

# Process data rows: copy value from column 78 to first column
for row in modified_rows:
    # Create copy of row
    new_row = row.copy()
    
    # Get connection string (column 78)
    connected_str = new_row[-1]
    
    # Replace first column value with connection string
    new_row[0] = connected_str
    
    # Add to final output list
    final_output_rows.append(new_row[:77])  # Keep only first 77 columns

# Step 6: Append processed data to output_file1
with open(output_file1, 'a', newline='') as f_out1_append:
    writer = csv.writer(f_out1_append)
    
    # Skip header row, only write data rows
    for row in final_output_rows[1:]:
        writer.writerow(row)

# Step 7: Update junction-refinement-result file (with replaced first column)
with open(output_file3, 'w', newline='') as f_out3:
    writer = csv.writer(f_out3)
    # Write new header (including Connected column)
    writer.writerow(new_header)
    
    # Write modified data (including all 78 columns)
    for row in modified_rows:
        # Create copy of row
        new_row = row.copy()
        
        # Replace first column with connection string
        new_row[0] = row[-1]
        
        writer.writerow(new_row)

print("Processing completed!")
print(f"Results saved to: {output_file1}, {output_file2}, {output_file3}")