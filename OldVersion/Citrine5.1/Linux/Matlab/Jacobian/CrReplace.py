import os
import shutil

input_folders = ['Crout(Fj1)', 'Crout(Fj0)']

try:
    for input_folder in input_folders:
        output_folder = f'Crout(Fj{os.path.basename(input_folder)[6:7]})'

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        for root, dirs, files in os.walk(input_folder):
            for filename in files:
                if filename.endswith('.txt'):
                    input_file_path = os.path.join(root, filename)
                    relative_path = os.path.relpath(input_file_path, input_folder)
                    output_file_path = os.path.join(output_folder, relative_path)

                    output_file_dir = os.path.dirname(output_file_path)
                    if not os.path.exists(output_file_dir):
                        os.makedirs(output_file_dir)

                    with open(input_file_path, 'r', encoding='utf-8') as file:
                        content = file.read()

                    variables = [
                        ('u0old', 'Yold(i*10+ 0)'), ('v0old', 'Yold(i*10+ 1)'), ('w0old', 'Yold(i*10+ 2)'), ('T0old', 'Yold(i*10+ 3)'),
                        ('Sn0old', 'Yold(i*10+ 4)'), ('Sb0old', 'Yold(i*10+ 5)'), ('Theta0old', 'Yold(i*10+ 6)'), ('Phi0old', 'Yold(i*10+ 7)'),
                        ('O2mega0old', 'Yold(i*10+ 8)'), ('O3mega0old', 'Yold(i*10+ 9)'), ('u1old', 'Yold(i*10+ 10)'), ('v1old', 'Yold(i*10+ 11)'),
                        ('w1old', 'Yold(i*10+ 12)'), ('T1old', 'Yold(i*10+ 13)'), ('Sn1old', 'Yold(i*10+ 14)'), ('Sb1old', 'Yold(i*10+ 15)'),
                        ('Theta1old', 'Yold(i*10+ 16)'), ('Phi1old', 'Yold(i*10+ 17)'), ('O2mega1old', 'Yold(i*10+ 18)'), ('O3mega1old', 'Yold(i*10+ 19)'),
                        ('u0new', 'Ynew(i*10+ 0)'), ('v0new', 'Ynew(i*10+ 1)'), ('w0new', 'Ynew(i*10+ 2)'), ('T0new', 'Ynew(i*10+ 3)'),
                        ('Sn0new', 'Ynew(i*10+ 4)'), ('Sb0new', 'Ynew(i*10+ 5)'), ('Theta0new', 'Ynew(i*10+ 6)'), ('Phi0new', 'Ynew(i*10+ 7)'),
                        ('O2mega0new', 'Ynew(i*10+ 8)'), ('O3mega0new', 'Ynew(i*10+ 9)'), ('u1new', 'Ynew(i*10+ 10)'), ('v1new', 'Ynew(i*10+ 11)'),
                        ('w1new', 'Ynew(i*10+ 12)'), ('T1new', 'Ynew(i*10+ 13)'), ('Sn1new', 'Ynew(i*10+ 14)'), ('Sb1new', 'Ynew(i*10+ 15)'),
                        ('Theta1new', 'Ynew(i*10+ 16)'), ('Phi1new', 'Ynew(i*10+ 17)'), ('O2mega1new', 'Ynew(i*10+ 18)'), ('O3mega1new', 'Ynew(i*10+ 19)')
                    ]

                    updated_content = content
                    for old, new in variables:
                        updated_content = updated_content.replace(old, new)

                    with open(output_file_path, 'w', encoding='utf-8') as file:
                        file.write(updated_content)

                    print(f"Processed {input_file_path} -> {output_file_path}")

        shutil.rmtree(input_folder)
        print(f"Deleted the original folder: {input_folder}")

    print("Done! All files in the output folders have been processed and the original folders have been deleted.")

except Exception as e:
    print(f"An error occurred: {e}")