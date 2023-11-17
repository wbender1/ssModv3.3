# Immunological Shape Space Computational Model for Respiratory Syncytial Viruses
ssModv3.3 Runs on Python 3.11

ssModv3.3.py is the main script to be executed.

config.py is a configuration file which is used to set variables for the simulation. 

The default simulation is configured to administer 3 vaccination doses at 4 weeks, 8 weeks, and 12 weeks. A challenge virus will also be administered 4 weeks after the last vaccination dose and the simulation will run for another 4 weeks.

The vaccination doses and challenge virus are chosen from a list of 12 RSV strains that are comprised of 4 G-protein Epitopes and 6 F-protein Epitopes. 

Each vaccination dose is split in two. This allows for two different strains to be administered together as a Bivalent or Heterologous vaccination. 
