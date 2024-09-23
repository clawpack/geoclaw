This directory contains a few examples to get users started using
this software.  In most directories you should be able to type:
    make all
and the code should run and a set of plots produced.  

However, see the individual README files in the directories.

To run all examples in all subdirectories:
    python $CLAW/clawutil/src/python/clawutil/run_examples.py
You might want to first check your environment variables, or you can set
them explicitly by modifying the script.

In addition to doing `make all`, this script also runs any Jupyter notebooks
found in the subdirectories.

It may take some time for these examples to all run. Two files 
    make_all_output.txt
    make_all_errors.txt
will be created that contain the output normally sent to the screen and any
error messages from running the examples.

A quicker set of tests can be run from the geoclaw/tests directory.  See the
README file in that directory.

