from jinja2 import Environment, FileSystemLoader

def getparameters(vardict):
    paramdict={}
    for key,val in vardict.items():
        if not val is None:
            paramdict[key]=val
    return paramdict

def generate_report(plotlist,tranplot, parameters,
                    table, filename, pertranscript,):
    '''
    docstring
    '''

    env = Environment(loader=FileSystemLoader('Translonpredictor'))
    template = env.get_template('report.html')
    report = template.render(
        plotlist=plotlist,
        tranplot=tranplot,
        table=table,
        filename=filename,
        tranplotlist=pertranscript,
        parameters=parameters
        )
    #to save the results
    with open(f"report_{filename}.html", "w") as fh:
        fh.write(report)

    return