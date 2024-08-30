import os

def get_files_path(folder_path) -> list[str]:
    results = []
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".fcs"):
            file_path = os.path.join(folder_path, file_name)
            results.append(file_path)
    return results    
    