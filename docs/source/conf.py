from Shiba import __version__

version = __version__

project = 'CHECK_PACKAGE_NAME'
copyright = f"rtuovine, DoryeLEsteras"
author = 'rtuovine, DoryeLEsteras'

extensions = []

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'renku'
html_static_path = ['_static']

variables_to_export = [
    "version",
]

html_context = {
    "default_mode": "light",
    "display_github": False,  
    "github_user": "DoryeLEsteras",  
    "github_repo": "Shiba",  
    "doc_path": "docs/source", 
}

