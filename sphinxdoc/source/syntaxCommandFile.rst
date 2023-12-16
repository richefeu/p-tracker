The command file 
================


In p-tracker, the commands and parameters are set in a command file, whilch is a simple text file.
The file contains specific keywords that define various aspects of the configuration. 
These keywords are used throughout the file to specify different settings and parameters.
Comments can added, line by line, when the caracter `#`, `!` or `/` is found.
This page provides a description of the keywords and syntax for the main usage of the tool: the particle image tracking. For this, the first usefull keywords is:

- ``procedure`` (*string*) **procedure_name**

  Set the procedure to be executed (default is `particle_tracking`).


Multi-threading
---------------

- ``wanted_num_threads`` (*int*) **value**

  Number of threads to be used for parallel computing (maximum is 48).

Images
------

- ``image_name`` (*string*) **c-style-path**

  Path (relative or not) of numbered image files. The number is set in a c-style printf formatName of the numbered images, e.g., `../ImageData/TEST%02d.tiff` 

Tracked point positions
-----------------------

- ``grains_to_follow`` (*string*) **filename**

  Todo.






