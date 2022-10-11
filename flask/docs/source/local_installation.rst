.. _local_installation:

####################################
Local installation of LocusFocus
####################################

LocusFocus may be cloned from the `GitHub repository <https://github.com/naim-panjwani/LocusFocus>`_ 
and run as a local webserver in Mac or Linux distributions. 

The GTEx database, and 1000 Genomes PLINK files are optional. GTEx may be built using the 
`initdb_GTExV7.py <https://github.com/naim-panjwani/LocusFocus/blob/master/misc/initdb_GTExV7.py>`_ 
and `initdb_GTExV8.py <https://github.com/naim-panjwani/LocusFocus/blob/master/misc/initdb_GTExV7.py>`_ scripts.
The script pushes `GTEx eQTL summary statistics <https://gtexportal.org/home/datasets>`_ 
into a NoSQL `MongoDB <https://www.mongodb.com/>`_ database.

The app has been tested in Linux and Mac environments.
For Windows, you may use the `Windows Subsystem for Linux (WSL) <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_
to emulate a Linux environment. The app has been tested using the Ubuntu WSL system.

All required programs and packages are listed in a `yml file <https://github.com/naim-panjwani/LocusFocus/blob/master/environment.yml>`_ and a 
`conda spec file <https://github.com/naim-panjwani/LocusFocus/blob/master/spec-file.txt>`_ are also available for direct explicit
installation of the exact conda virtual environment the app requires for a successful run.

To clone and install the required packages, simply run:

.. code-block:: console
   :caption: *Clone and install required conda environment*

   > git clone https://github.com/naim-panjwani/LocusFocus.git
   > conda create --name <env> --file spec-file.txt

Next, to run the application locally, simply run the app.py script:

.. code-block:: console
   :caption: *Run the app*

   > python app.py

The above command will issue a web address to access the application locally. 
Simply type the address given in your browser.

