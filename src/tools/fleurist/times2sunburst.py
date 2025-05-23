import json
import plotly.express as px

def generate_sunburst_plot(json_file, output_file="juDFT_times_plot.html"):
    """
    Generate a sunburst plot from a JSON file.

    Parameters:
        json_file (str): Path to the JSON file containing hierarchical data.
        output_file (str): Path to save the generated sunburst plot as an HTML file.
    """
    # Load the JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)

    # Convert JSON data to a format suitable for Plotly
    def flatten_json(data, parent="",runtime=0):
        rows = []
        timer=data['timername']
        total=float(data['totaltime'])
        if 'ncalls' in data:
            ncalls=data['ncalls']
            variance=round((float(data['maxtime'])-float(data['mintime']))/float(data['maxtime'])*100,2)
        else:
            variance=0
            ncalls=1    
        if runtime==0: runtime=total
        percent=total/runtime*100
        current_path = f"{parent}/{timer}" if parent else timer
        if total>runtime*0.00001:
            rows.append({"id": current_path, "label": timer, "parent": parent, "runtime": total,"percent": f'{percent:.2f}%',"calls":ncalls, "variance": variance})
            if 'subtimers' in data:
                for subtimer in data['subtimers']:
                    rows.extend(flatten_json(subtimer, timer,runtime))
        
        return rows

    flat_data = flatten_json(data)

    # Create a DataFrame for the sunburst plot
    import pandas as pd
    df = pd.DataFrame(flat_data)
    
    # Generate the sunburst plot
    fig = px.sunburst(
        df,
        names="label",
        parents="parent",
        values="runtime",
        hover_data=["percent", "calls","variance"], color="variance",
        title="FLEUR runtime breakdown" ,branchvalues="total" )

    # Save the plot to an HTML file
    fig.write_html(output_file)
    print(f"Sunburst plot saved to {output_file}")
