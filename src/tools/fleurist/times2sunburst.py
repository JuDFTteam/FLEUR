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
    def flatten_json(data, parent="",runtime=0,sector_time=0,scaling_data=None):
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
        if runtime==0: 
            runtime=total
            sector_time=total
        if scaling_data:
            scaling=scaling_data['totaltime']/total
        percent=total/runtime*100
        current_path = f"{parent}/{timer}" if parent else timer
        total=min(total,sector_time)
        if total>runtime*0.01 and total>sector_time*0.01    :
            if timer == "Total Run":
                timer="" #reset parent timer
            else:
                rows.append({"id": current_path, "label": timer, "parent": parent, "runtime": total,"percent": f'{percent:.2f}%',"calls":ncalls, "variance": variance, "scaling": scaling})
                print(rows[-1])
            if 'subtimers' in data:
                if scaling_data:
                    sbtimer_scaling=scaling_data['subtimers']
                for subtimer in data['subtimers']:
                    if scaling_data:
                        subtimer_scale=sbtimer_scaling.pop(0)
                    else:
                        subtimer_scale=None  
                    rows.extend(flatten_json(subtimer, timer,runtime,total,subtimer_scale))
        
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
        values="runtime",color=color,
        hover_data=["percent", "calls","variance"], 
        title="FLEUR runtime breakdown" ,branchvalues="total" )
    fig.show()
    # Save the plot to an HTML file
    fig.write_html(output_file)
    print(f"Sunburst plot saved to {output_file}")
