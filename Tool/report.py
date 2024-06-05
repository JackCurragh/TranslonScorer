from jinja2 import Environment, FileSystemLoader

def generate_report(plotlist,tranplot,table):
    '''
    docstring
    '''

    env = Environment(loader=FileSystemLoader('Tool'))
    template = env.get_template('report.html')
    report = template.render(
        plotlist=plotlist,
        tranplot=tranplot,
        table=table
        )
    #to save the results
    with open("Report.html", "w") as fh:
        fh.write(report)

    return