# import os

# import jinja2
# from jinja2 import BaseLoader

# def generate_prediction_report(
#     output_file,
#     **template_data
# ):


#     with open(os.path.abspath(os.path.join(os.path.dirname(__file__), 'include/templates', 'report.html'))) as f:
#         template_str = f.read()

#     template = jinja2.Environment(loader=BaseLoader).from_string(template_str)
#     html = template.render(
#         **template_data
#     )

#     # Write the HTML file
#     with open(output_file, 'w') as f:
#         f.write(html)
