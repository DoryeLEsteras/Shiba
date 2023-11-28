from Shiba import __release_date__, __version__

version = __version__
release_date = __release_date__

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
    "release_date",
]

extlinks = {
    "githubDLE": (
        "https://github.com/adrybakov/rad-tools/tree/stable/docs/examples/%s",
        "%s",
    ),

    "numpy": (
        "https://numpy.org/doc/stable/reference/generated/numpy.%s.html",
        "numpy.%s",
    )
}

html_context = {
    "default_mode": "light",
    "display_github": False,  
    "github_user": "DoryeLEsteras",  
    "github_repo": "Shiba",  
    "doc_path": "docs/source", 
}

html_theme_options = {
    "collapse_navigation": True,
    "use_edit_page_button": True,
    "navbar_center": ["version-switcher", "navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "switcher": {
        "version_match": switcher_version,
        "json_url": "https://rad-tools.org/en/stable/_static/versions.json",
    },
    "navbar_align": "left",
    "logo": {
        "image_light": "_static/logo_black.png",
        "image_dark": "_static/logo_white.png",
    },
    "header_links_before_dropdown": 4,
    "icon_links": [
        {
            "name": "Twitter",
            "url": "https://twitter.com/adrybakov",
            "icon": "fa-brands fa-twitter",
        },
        {
            "name": "GitHub",
            "url": "https://github.com/adrybakov/rad-tools",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/rad-tools/",
            "icon": "fa-solid fa-box",
        },
    ],
}