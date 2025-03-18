import os
import re

def rename_files_arrange():
    # Get user input
    static_part = input("Enter the part of the name you want to keep: ")
    dynamic_part = input("If you want to arrange the files, type (arrange): ")
    file_extension = input("Enter the file extension to change (e.g., .gro, .txt), or leave empty for all: ")

    # Get a list of files in the directory
    files = [f for f in os.listdir() if os.path.isfile(f)]

    # Filter files by extension, if provided
    if file_extension:
        files = [f for f in files if f.endswith(file_extension)]

    # Filter files matching the given pattern
    matching_files = [f for f in files if f.startswith(static_part)]

    # Sort files based on the number in the dynamic part (if it exists)
    matching_files.sort(key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else 0)

    # Rename files in ascending order
    for i, filename in enumerate(matching_files, start=1):
        new_filename = f"{static_part}{i}{file_extension}"
        os.rename(filename, new_filename)
        print(f"Renamed: {filename} -> {new_filename}")

    print("Operation completed!")

if __name__ == "__main__":
    rename_files_arrange()

