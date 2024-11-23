import os
import glob
import re

def find_matching_parenthesis(expr, pos, direction):
    stack = 1
    if direction == 'left':
        for i in range(pos - 1, -1, -1):
            if expr[i] == ')':
                stack += 1
            elif expr[i] == '(':
                stack -= 1
                if stack == 0:
                    return i
    elif direction == 'right':
        for i in range(pos + 1, len(expr)):
            if expr[i] == '(':
                stack += 1
            elif expr[i] == ')':
                stack -= 1
                if stack == 0:
                    return i
    return -1

def replace_powers(expr):
    while '^' in expr:
        caret_pos = expr.find('^')
        
        if caret_pos == 0:
            raise ValueError("The expression starts with a caret (^), but the base number cannot be found.")
        
        if expr[caret_pos - 1] == ')':
            base_end = caret_pos - 1
            base_start = find_matching_parenthesis(expr, base_end, 'left')
            if base_start == -1:
                raise ValueError("The parentheses do not match, and the starting position of the base cannot be found.")
            base = expr[base_start:base_end + 1]
        else:
            base_end = caret_pos - 1
            base_start = base_end
            while base_start >= 0 and expr[base_start] not in '+-*/%() ':
                base_start -= 1
            base_start += 1
            base = expr[base_start:base_end + 1]
        
        exponent_start = caret_pos + 1
        if exponent_start >= len(expr):
            raise ValueError("The ^ symbol is followed by a missing exponent.")
        
        if expr[exponent_start] == '(':
            exponent_end = find_matching_parenthesis(expr, exponent_start, 'right')
            if exponent_end == -1:
                raise ValueError("Mismatched parentheses, unable to find the end position of the exponent.")
            exponent = expr[exponent_start:exponent_end + 1]
        else:
            exponent_end = exponent_start
            while exponent_end < len(expr) and expr[exponent_end] not in '+-*/%() ':
                exponent_end += 1
            exponent = expr[exponent_start:exponent_end]
        
        # Check if exponent is 0.5 or 1/2
        if exponent == '0.5' or exponent == '(1/2)' or exponent == '1/2':
            result_str = f"sqrt({base})"
        else:
            result_str = f"pow({base}, {exponent})"
        
        if expr[exponent_start] == '(':
            expr = expr[:base_start] + result_str + expr[exponent_end + 1:]
        else:
            expr = expr[:base_start] + result_str + expr[exponent_end:]
    
    return expr

def process_file(file_path, output_dir=None):
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    new_lines = []
    changes_made = False
    for line in lines:
        if '^' in line:
            try:
                new_line = replace_powers(line)
                if new_line != line:
                    changes_made = True
                new_lines.append(new_line)
            except ValueError as e:
                print(f"An error occurred while processing lines in file {file_path}: {e}")
                new_lines.append(line)
        else:
            new_lines.append(line)
    
    if output_dir:
        base_name = os.path.basename(file_path)
        new_file_path = os.path.join(output_dir, base_name)
        with open(new_file_path, 'w', encoding='utf-8') as file:
            file.writelines(new_lines)
        if changes_made:
            print(f"A new file has been created: {new_file_path}")
    else:
        if changes_made:
            with open(file_path, 'w', encoding='utf-8') as file:
                file.writelines(new_lines)
            print(f"File modified: {file_path}")
        else:
            print(f"File not modified: {file_path}")


def main():
    input_folders = ['FJ1 output_txt_files', 'FJ0 output_txt_files']
    output_folders = ['Crout(Fj1)', 'Crout(Fj0)']

    for input_dir, output_dir in zip(input_folders, output_folders):
        if not os.path.exists(input_dir):
            print(f"Directory not found: {input_dir}")
            continue

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        txt_files = glob.glob(os.path.join(input_dir, '*.txt'))
        
        if not txt_files:
            print(f"No .txt files found in the directory {input_dir}.")
            continue

        for txt_file in txt_files:
            print(f"Processing file: {txt_file}")
            process_file(txt_file, output_dir)
    
    print("All files have been processed.")

if __name__ == "__main__":
    main()