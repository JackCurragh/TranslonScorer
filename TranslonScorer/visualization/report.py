"""
Report generation functionality for TranslonScorer.

This module handles the generation of HTML reports from analysis results,
using Jinja2 templates to create interactive visualizations.
"""

from jinja2 import Environment, FileSystemLoader
from pathlib import Path
from ..utils.logging import log_info, log_warning


def getparameters(vardict):
    """
    Filter out None values from a dictionary.

    Args:
        vardict (dict): Dictionary with parameter names and values

    Returns:
        dict: Filtered dictionary containing only non-None values
    """
    paramdict = {}
    for key, val in vardict.items():
        if val is not None:
            paramdict[key] = val
    return paramdict


def generate_report(plotlist, tranplot, parameters, table, filename, pertranscript):
    """
    Generate an HTML report using a template.

    Args:
        plotlist (list): List of HTML strings containing metagene plots
        tranplot (list): List of HTML strings containing transcript plots
        parameters (dict): Dictionary of analysis parameters
        table (str): HTML string containing summary table
        filename (str): Output filename for report
        pertranscript (list): List of HTML strings with transcript-specific plots
    """
    log_info(f"Generating report {filename}_report.html...")

    # Get package directory for template loading
    pkg_dir = Path(__file__).parent.parent
    template_dir = pkg_dir / "templates"
    template_dir.mkdir(exist_ok=True)

    # Create template if it doesn't exist
    template_path = template_dir / "report.html"
    if not template_path.exists():
        create_report_template(template_path)

    # Generate report
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template("report.html")

    # Filter parameters
    if parameters:
        parameters = getparameters(parameters)

    report = template.render(
        plotlist=plotlist,
        tranplot=tranplot,
        table=table,
        filename=filename,
        tranplotlist=pertranscript,
        parameters=parameters,
    )

    # Save report
    output_path = Path(filename).parent
    output_path.mkdir(parents=True, exist_ok=True)

    with open(f"{filename}_report.html", "w") as fh:
        fh.write(report)

    log_info("Report generated successfully")


def create_report_template(template_path):
    """
    Create the default report template.

    Args:
        template_path (Path): Path where template should be created
    """
    template = """<!DOCTYPE html>
<html>
<head>
    <title>TranslonScorer Report - {{filename}}</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 40px;
            line-height: 1.6;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
        }
        .section {
            margin-bottom: 40px;
            padding: 20px;
            background: #f9f9f9;
            border-radius: 5px;
        }
        h1, h2 {
            color: #333;
        }
        .parameters {
            background: #eef;
            padding: 15px;
            border-radius: 5px;
        }
        .plot {
            margin: 20px 0;
            padding: 10px;
            background: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>TranslonScorer Analysis Report</h1>
        <div class="section">
            <h2>Analysis Parameters</h2>
            <div class="parameters">
                {% if parameters %}
                    {% for key, value in parameters.items() %}
                        <p><strong>{{ key }}:</strong> {{ value }}</p>
                    {% endfor %}
                {% else %}
                    <p>No parameters provided</p>
                {% endif %}
            </div>
        </div>

        <div class="section">
            <h2>Summary Table</h2>
            {{ table }}
        </div>

        <div class="section">
            <h2>Metagene Profiles</h2>
            {% for plot in plotlist %}
                <div class="plot">{{ plot }}</div>
            {% endfor %}
        </div>

        <div class="section">
            <h2>Transcript Summary</h2>
            <div class="plot">{{ tranplot }}</div>
        </div>

        <div class="section">
            <h2>Individual Transcript Plots</h2>
            {% for plot in tranplotlist %}
                <div class="plot">{{ plot }}</div>
            {% endfor %}
        </div>
    </div>
</body>
</html>"""

    with open(template_path, "w") as f:
        f.write(template)
