import os

# replace with own paths to executables if not in PATH

TABULATE_ALIGNMENTS_EXECUTABLE = os.environ.get(
    "TABULATE_ALIGNMENTS_EXECUTABLE", "tabulate-alignments") 
CLIQUESNV_JAR = os.environ.get("CLIQUESNV_JAR", "clique-snv.jar")
REMOVECLIPPING_SCRIPT = os.environ.get("REMOVECLIPPING_SCRIPT", "removeclipping_py3.py")
