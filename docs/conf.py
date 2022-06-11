"""Sphinx configuration."""
from datetime import datetime


project = "CsfDock"
author = "Rahul Brahma"
copyright = f"{datetime.now().year}, {author}"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_click",
    "sphinx_rtd_theme",
]
autodoc_typehints = "description"
html_theme = "sphinx_rtd_theme"


import os

extensions = [... "rtds_action"]

# The name of your GitHub repository
rtds_action_github_repo = "USERNAME/REPONAME"

# The path where the artifact should be extracted
# Note: this is relative to the conf.py file!
rtds_action_path = "tutorials"

# The "prefix" used in the `upload-artifact` step of the action
rtds_action_artifact_prefix = "notebooks-for-"

# A GitHub personal access token is required, more info below
rtds_action_github_token = os.environ["GITHUB_TOKEN"]
