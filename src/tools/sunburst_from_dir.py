import os
import plotly.express as px
import pandas as pd

def get_directory_structure(path):
    """
    Recursively get the directory structure and count the number of lines in files.

    Parameters:
        path (str): The root directory path.

    Returns:
        list: A list of dictionaries containing 'id', 'parent', and 'lines'.
    """
    data = []
    valid_extensions = {".f90", ".f", ".F90", ".F"}

    def count_lines(file_path):
        """Count the number of lines in a file."""
        try:
            with open(file_path, 'r') as file:
                return sum(1 for _ in file)
        except Exception:
            return 0

    def traverse_directory(current_path, parent=""):
        total_lines = 0
        for entry in os.scandir(current_path):
            entry_path = os.path.join(current_path, entry.name)
            if entry.is_dir():
                dir_lines = traverse_directory(entry_path, entry.name)
                data.append({"id": entry.name, "parent": parent, "lines": dir_lines})
                total_lines += dir_lines
            elif os.path.splitext(entry.name)[1] in valid_extensions:
                file_lines = count_lines(entry_path)
                data.append({"id": entry.name, "parent": parent, "lines": file_lines})
                total_lines += file_lines
        return total_lines

    root_lines = traverse_directory(path,path)
    data.append({"id": path, "parent": "", "lines": root_lines})
    return data

def plot_sunburst(directory, output_file="directory_sunburst.html"):
    """
    Generate a sunburst plot for the directory structure based on the number of lines in files.

    Parameters:
        directory (str): The root directory to analyze.
        output_file (str): The output HTML file for the sunburst plot.
    """
    # Get the directory structure
    structure = get_directory_structure(directory)

    # Convert to a DataFrame
    df = pd.DataFrame(structure)

    # Generate the sunburst plot
    fig = px.sunburst(
        df,
        names="id",
        parents="parent",
        values="lines",branchvalues="total",
        title=f"Sunburst Diagram of {directory} (by number of lines)",
    )

    # Show the plot 
    fig.show()
    #fig.write_html(output_file)
    #print(f"Sunburst plot saved to {output_file}")

# Example usage
if __name__ == "__main__":
    directory_to_analyze = "."  # Change this to the directory you want to analyze
    plot_sunburst(directory_to_analyze)