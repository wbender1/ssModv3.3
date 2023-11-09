#Configurations for ssModv3 Simulation
#Doses at 4 weeks, 8 weeks, 12 weeks (28d, 56d, 84d)
vacc_interval1 = 28.0
vacc_interval2 = 28.0
vacc_interval3 = 28.0
#1 dose of Challenge Virus Infection at week 16 (112d)
infect_interval = 28.0
infect = "yes"
challengeVirus = 19
#Simulation ends after week 20 (140d)
stop_interval = 28.0
#Initial Antigen amounts for each dose and infection
ag1_initialA = 180
ag1_initialB = 180
ag2_initialA = 180
ag2_initialB = 180
ag3_initialA = 180
ag3_initialB = 180
virus_initial = 180
#Antigen numbers for each dose
prime_agA = 7
prime_agB = 7
boost1_agA = 7
boost1_agB = 7
boost2_agA = 7
boost2_agB = 7
#Neutralization configuration
AG1_ep1_neutralization = 1.0
AG1_ep2_neutralization = 1.0
AG1_ep3_neutralization = 1.0
AG1_ep4_neutralization = 1.0
AG2_ep1_neutralization = 1.0
AG2_ep2_neutralization = 1.0
AG2_ep3_neutralization = 1.0
AG2_ep4_neutralization = 1.0
AG2_ep5_neutralization = 1.0
AG2_ep6_neutralization = 1.0
#Immunogenicity configuration
#(a value of 1.0 for AG1 will cause the simulation to use G only Antigens)
#(a value of 1.0 for AG2 will cause the simulation to use F only Antigens)
AG1_ep1_immunogenicity = 1.0
AG1_ep2_immunogenicity = 1.0
AG1_ep3_immunogenicity = 1.0
AG1_ep4_immunogenicity = 1.0
AG2_ep1_immunogenicity = 0.0
AG2_ep2_immunogenicity = 0.0
AG2_ep3_immunogenicity = 0.0
AG2_ep4_immunogenicity = 0.0
AG2_ep5_immunogenicity = 0.0
AG2_ep6_immunogenicity = 0.0
#Antigen Clearance configuration
AG1_ep1_clearance = 1.0
AG1_ep2_clearance = 1.0
AG1_ep3_clearance = 1.0
AG1_ep4_clearance = 1.0
AG2_ep1_clearance = 1.0
AG2_ep2_clearance = 1.0
AG2_ep3_clearance = 1.0
AG2_ep4_clearance = 1.0
AG2_ep5_clearance = 1.0
AG2_ep6_clearance = 1.0
#Other configurations
ag_decay = 12.0
stimulation = 1.0
clearance = 0.0001
Bcell_initial = 50000000
Bcell_carry = 5000
Bcell_aff = 10.0
AB_aff = 2.5
tau = 8.0
