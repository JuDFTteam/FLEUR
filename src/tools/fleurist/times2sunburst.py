import json
import plotly.express as px

def generate_sunburst_plot(json_file, output_file="juDFT_times_plot.html",scalingFile=None):
    """
    Generate a sunburst plot from a JSON file.

    Parameters:
        json_file (str): Path to the JSON file containing hierarchical data.
        output_file (str): Path to save the generated sunburst plot as an HTML file.
    """
    # Load the JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)

    if scalingFile:
        # Load the scaling data
        with open(scalingFile, 'r') as file:
            scaling_data = json.load(file)
    else:
        scaling_data= None

    # Convert JSON data to a format suitable for Plotly
    def flatten_json(data, parent="",runtime=0,scaling_data=None):
        scaling=0
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
        if scaling_data:
            scaling=scaling_data['totaltime']/total
        percent=total/runtime*100
        current_path = f"{parent}/{timer}" if parent else timer
        if total>runtime*0.001:
            rows.append({"id": current_path, "label": timer, "parent": parent, "runtime": total,"percent": f'{percent:.2f}%',"calls":ncalls, "variance": variance, "scaling": scaling})
            if 'subtimers' in data:
                if scaling_data:
                    sbtimer_scaling=scaling_data['subtimers']
                for subtimer in data['subtimers']:
                    if scaling_data:
                        subtimer_scale=sbtimer_scaling.pop(0)
                    else:
                        subtimer_scale=None  
                    rows.extend(flatten_json(subtimer, timer,runtime,subtimer_scale))
        
        return rows

    flat_data = flatten_json(data,scaling_data=scaling_data)

    # Create a DataFrame for the sunburst plot
    import pandas as pd
    df = pd.DataFrame(flat_data)
    
    color="scaling" if scaling_data else "variance"

    # Generate the sunburst plot
    fig = px.sunburst(
        df,
        names="label",
        parents="parent",
        values="runtime",
        hover_data=["percent", "calls","variance"], color=color,
        title="FLEUR runtime breakdown" ,branchvalues="total" )

    # Save the plot to an HTML file
    fig.write_html(output_file)
    print(f"Sunburst plot saved to {output_file}")
