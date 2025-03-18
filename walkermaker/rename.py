import os

def rename_files():
    # Get user input
    original_text = input("Enter the part of the name to replace: ")
    replacement_text = input("Enter the new text (or leave empty to remove): ")
    file_extension = input("Enter the file extension to modify (e.g., .gro, .txt), or leave empty for all: ")

    # Get a list of files in the current directory
    files = [f for f in os.listdir() if os.path.isfile(f)]

    # Filter files by extension, if provided
    if file_extension:
        files = [f for f in files if f.endswith(file_extension)]

    # Iterate through files and rename them
    for filename in files:
        if original_text in filename:
            new_filename = filename.replace(original_text, replacement_text)
            os.rename(filename, new_filename)
            print(f"Renamed: {filename} -> {new_filename}")

    print("Operation completed!")

if __name__ == "__main__":
    rename_files()

