import click
import sys
import fleurjobs
import fleuristconf
import fleurinput

@click.group()
def cmd_main():
    """
    FLEURist (not Florist) is a tool to generate parallelization suggestions for the FLEUR code. It should be called 
    from the directory in which the FLEUR inp.xml file resides. This file will be evaluated and suggestions for parallelization
    generated. The tool can also directly submit to a slurm queue or run the code on the local computer.

    Before using FLEURist you should configure the computers as described by 'FLEURist computer add --help'.
    """
    pass

@cmd_main.command('suggest')
@click.option("-c","--computer",help="Choose name of the computer/queue to use",default=None)
@click.option("-n","--nodes",help="Choose number of nodes to use",default=0)
@click.option("-e","--efficiency",help="Minimal parallel efficiency expected for multi-node jobs",default=0.5)
def suggest(computer,nodes,efficiency):
    """
    Print suggestions for possible parallel setups for FLEUR.
    """
    try:  
        input=fleurinput.read_inpxml("./inp.xml")
    except:
        click.secho("Could not read inp.xml file",fg='red')
        click.secho("Please check that the file exists and is a valid FLEUR input-file.")    
        import sys
        sys.exit()
    click.secho("Suggesting job settings for inp.xml in current directory",fg="red")
    click.secho("========================================================")
    click.secho(f"Number of k-points: {input['nkpt']}")
    click.secho(f"Size of the setup : {input['basis']}")
    click.secho("========================================================\n")
    machine,resources=fleurjobs.suggest(input,computer,int(nodes),float(efficiency))
    if machine:
        for i in range(len(machine)):
            if i==0:
                click.secho(f"0    Default suggestion:{machine[i].name}",fg='red')
                click.secho("-----===================",fg='red')
            else:
                click.secho(f"{i}.  Alternative suggestion:{machine[i].name}")
                click.secho("-----========================")
            click.secho(resources[i].description())
    else:
        click.secho("No suggestion found.",fg='red')
        if nodes:
            click.secho("Try to relax the required number of nodes")
        else:
            click.secho("Please report if you feel that your setup should be ok for FLEUR")    


@cmd_main.command('submit')
@click.option("-c","--computer",help="Choose name of the computer/queue",default=None)
@click.option("-n","--nodes",help="Choose number of nodes to use",default=0)
@click.option("-e","--efficiency",help="Minimal parallel efficiency expected for multi-node jobs",default=0.5)
@click.option("--create",is_flag=True,help="Create a jobfile, but do not submit it",default=False)
@click.option("-s","--suggestion",help="Number of the suggested config to use (see the suggest command)",default=0)
def submit(computer,nodes,efficiency,create,suggestion):
    """
    Submit a job. Either to slurm or run it directly on the localhost. It can be advisable to use 'suggest' before to check the suggested config.
    """  
    try:  
        input=fleurinput.read_inpxml("./inp.xml")
    except:
        click.secho("Could not read inp.xml file",fg='red')
        click.secho("Please check that the file exists and is a valid FLEUR input-file.")    
        sys.exit()

    machine,resources=fleurjobs.suggest(input,computer,nodes,efficiency)

    if machine:
        if suggestion+1>len(machine) or suggestion<0:
            click.secho("You have choosen a suggestion number not available.", fg='red')
            click.secho("Please check 'FLEURist suggest' again for valid choices")
            sys.exit()           
        try:
            machine[suggestion].submit(resources[suggestion],create)
        except:
            click.secho("Job submission failed",fg='red')
            sys.exit()    
    else:
        click.secho("NO possible config found. Try different settings and the 'suggest' command before",fg='red')


@cmd_main.group('computer')
def cmd_config():
    """
    The 'Computer' command group contains commands to manage the computer config of FLEURist
    """

@cmd_config.command('add')
def add():
    """
    Add a computer to your configuration, at least a single computer is needed to be able to use FLEURist.
    This command will start an interactive process to determine the configuration of the computer.
    """  
    fleuristconf.add()

@cmd_config.command('remove')
@click.option("-n","--name",help="Name of the computer to remove",prompt="Name:")
def remove(name):
    """
    Remove a computer from your configuration.
    """  
    fleuristconf.remove(name)

@cmd_config.command('list')
def list():
    """
    List all computers known to your FLEURist install.
    """  
    fleuristconf.list()









if __name__ == '__main__':
    cmd_main()
