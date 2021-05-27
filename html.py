"""
Functions that aid in the production of the HTML webpages.
"""

from jinja2 import Environment, PackageLoader, FileSystemLoader, select_autoescape
from time import strftime
import unyt
from swiftsimio import load

class PlotsInPipeline:
    def __init__(self):
        self.plots_details = []

    def load_plots(self, title, caption, filename, id):
        some_details = dict(title=title, caption=caption, filename=filename, hash=id,)
        self.plots_details.append(some_details)

    def reset_plots_list(self):
        self.plots_details = []

def add_metadata_to_web(web, siminfo):
    data = load(siminfo.snapshot)

    TEMPLATE_FILE = "description.html"
    description_template = web.environment.get_template(TEMPLATE_FILE)
    web.variables["runs"].append(dict(description=description_template.render(data=data), ))
    return

def make_web(siminfo):
    data = load(siminfo.snapshot)
    web = WebpageCreator()
    web.add_run_metadata(data=data)
    return web

def add_web_section(web,title,caption,id,plots):
    web.variables["sections"].append(dict(title=title, caption=caption, id=id, plots=plots,))
    return

def render_web(web,output_path):
    web.add_metadata(page_name="SIDMPlots Page")
    web.render_webpage()
    web.save_html(f"{output_path}/index.html")
    return

def render_population_web(web,output_path):
    web.add_metadata(page_name="SIDMPlots Page")
    web.render_webpage()
    web.save_html(f"{output_path}/population.html")
    return

def format_number(number):
    """
    Formats a number from float (with or without units) to a latex-like number.
    """

    try:
        units = f"\\; {number.units.latex_repr}"
    except:
        units = ""

    try:
        mantissa, exponent = ("%.3g" % number).split("e+")
        exponent = f" \\times 10^{{{int(exponent)}}}"
    except:
        mantissa = "%.3g" % number
        exponent = ""

    return f"\\({mantissa}{exponent}{units}\\)"


def get_if_present_float(dictionary, value: str, input_unit=None, output_unit=None):
    """
    A replacement for .get() that also formats the number if present.
    
    Assumes data should be a float.
    """

    try:
        value = float(dictionary[value])

        if input_unit is not None:
            value = unyt.unyt_quantity(value, input_unit)

            if output_unit is not None:
                value.convert_to_units(output_unit)

        return format_number(value)
    except KeyError:
        return ""


def get_if_present_int(dictionary, value: str, input_unit=None, output_unit=None):
    """
    A replacement for .get() that also formats the number if present.
    
    Assumes data should be an integer.
    """

    try:
        value = int(dictionary[value])

        if input_unit is not None:
            value = unyt.unyt_quantity(value, input_unit)

            if output_unit is not None:
                value.convert_to_units(output_unit)

        return format_number(value)
    except KeyError:
        return ""


class WebpageCreator(object):
    """
    Creates webpages based on the information that is provided in
    the plots metadata through the autoplotter and the additional
    plotting interface provided through the pipeline.
    """

    variables: dict
    html: str

    def __init__(self):
        """
        Sets up the ``jinja`` templating system.
        """
        
        self.loader = FileSystemLoader(searchpath="./templates/")
        self.environment = Environment(loader=self.loader,autoescape=select_autoescape(["js"]))

        # Initialise empty variables dictionary, with the versions of
        # this package and the velociraptor package used.
        self.variables = dict(
            #pipeline_version=[],
            #velociraptor_version=[],
            creation_date=strftime(r"%Y-%m-%d"),
            sections=[],
            runs=[],
        )

        return

    def render_webpage(self):
        """
        Renders a webpage based on the internal variables stored in
        the ``variables`` dictionary.

        Parameters
        ----------

        template: str
            The name of the template that you wish to use. Defaults to
            "plot_viewer.html".

        Returns
        -------

        html: str
            The resulting HTML. This is also stored in ``.html``.
        """
        TEMPLATE_FILE = "plot_viewer.html"

        template = self.environment.get_template(TEMPLATE_FILE, parent="base.html")
        self.html = template.render(**self.variables)

        return self.html
    
    def add_run_metadata(self, data):
        """
            Adds the "run" metadata (using the user-defined description.html).
            
            Parameters
            ----------
            
            """
        
        TEMPLATE_FILE = "description.html"

        self.environment.filters["format_number"] = format_number
        self.environment.filters["get_if_present_float"] = get_if_present_float
        self.environment.filters["get_if_present_int"] = get_if_present_int

        description_template = self.environment.get_template(TEMPLATE_FILE)
        self.variables["runs"] = [dict(description=description_template.render(data=data),)]
        
        return


    def add_metadata(self, page_name: str):
        """
        Add additional metadata to the page.

        Parameters
        ----------

        page_name: str
            Name to put in the page title.
        """

        self.variables.update(dict(page_name=page_name))


    def save_html(self, filename: str):
        """
        Saves the html in ``self.html`` to the filename provided.

        Parameters
        ----------

        filename: str
            Full filename (including file path) to save the HTML as.
        """

        with open(filename, "w") as handle:
            handle.write(self.html)
