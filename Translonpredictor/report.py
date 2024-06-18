from jinja2 import Environment, FileSystemLoader


def getparameters(vardict):
    """
    Filter out None values from a dictionary and return a new dictionary with only non-None values.

    Parameters:
    - vardict (dict): A dictionary where keys are parameter names and values are parameter values.

    Returns:
    - paramdict (dict): A filtered dictionary containing only key-value pairs from `vardict` where the value is not None.

    This function iterates through each key-value pair in `vardict`. If the value associated with a key is not None, it
    adds that key-value pair to a new dictionary `paramdict`. Finally, it returns `paramdict` containing only non-None
    values from `vardict`.

    Example:
    >>> vardict = {'param1': 100, 'param2': None, 'param3': 'value'}
    >>> getparameters(vardict)
    {'param1': 100, 'param3': 'value'}

    Note: This function is useful for extracting valid parameters from a larger dictionary, especially when dealing with
    optional parameters that may have default None values.
    """
    paramdict = {}
    for key, val in vardict.items():
        if not val is None:
            paramdict[key] = val
    return paramdict


def generate_report(
    plotlist,
    tranplot,
    parameters,
    table,
    filename,
    pertranscript,
):
    """
    Generate an HTML report using a template, incorporating various plots, tables, and parameters.

    Parameters:
    - plotlist (list): List of HTML strings containing metagene plots and other plots.
    - tranplot (list): List of HTML strings containing individual transcript plots.
    - parameters (dict): Dictionary containing parameters and settings used in the analysis.
    - table (str): HTML string containing a table summarizing data from the analysis.
    - filename (str): Name of the output file for the generated report.
    - pertranscript (list): List of HTML strings containing plots specific to each transcript.

    Returns:
    - None

    This function takes several components needed for the report generation:
    - `plotlist`: Contains metagene plots and other plots for different types of ORFs.
    - `tranplot`: Contains individual transcript plots for top 10 ORFs per type.
    - `parameters`: Provides parameters and settings used in the analysis.
    - `table`: Represents a table summarizing data from the analysis.
    - `pertranscript`: Includes plots specific to each transcript.

    Using the Jinja2 templating engine, it loads an HTML template (`report.html`) from the 'Translonpredictor' directory.
    The function then renders the template with the provided components:
    - `plotlist`, `tranplot`, `table`, `filename`, `pertranscript`, and `parameters`.

    After rendering, it saves the generated HTML report as `report_<filename>.html` in the current working directory.

    Note: Ensure the 'Translonpredictor' directory contains the necessary `report.html` template for this function to work.

    Example usage:
    >>> plotlist = ['<div>Metagene plot 1</div>', '<div>Metagene plot 2</div>']
    >>> tranplot = ['<div>Transcript plot 1</div>', '<div>Transcript plot 2</div>']
    >>> parameters = {'param1': 100, 'param2': 'value'}
    >>> table = '<table><tr><td>...</td></tr></table>'
    >>> pertranscript = ['<div>Transcript-specific plot 1</div>', '<div>Transcript-specific plot 2</div>']
    >>> generate_report(plotlist, tranplot, parameters, table, 'output', pertranscript)

    This function assumes familiarity with the Jinja2 templating syntax and requires the 'report.html' template to be
    correctly structured to accept the provided components.
    """
    env = Environment(loader=FileSystemLoader("Translonpredictor"))
    template = env.get_template("report.html")
    report = template.render(
        plotlist=plotlist,
        tranplot=tranplot,
        table=table,
        filename=filename,
        tranplotlist=pertranscript,
        parameters=parameters,
    )
    # to save the results
    with open(f"report_{filename}.html", "w") as fh:
        fh.write(report)
    return
